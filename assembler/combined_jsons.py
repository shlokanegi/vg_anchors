#!/usr/bin/env python3
"""
Combine JSON files per chromosome based on chunk information from TSV file.

This script:
1. Reads a TSV file with chunk information (start_node, end_node, path_name)
2. Extracts chromosome names from path_name column
3. Groups chunks by chromosome
4. Finds and combines JSON files for chunks of each chromosome
"""

import json
import re
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Any, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing


def extract_chr_name(path_name: str) -> str:
    """
    Extract chromosome name from path_name.
    
    Examples:
        "CHM13#0#chr2" -> "chr2"
        "recombination#16#chr4#1" -> "chr4"
        "CHM13#0#chrX" -> "chrX"
        "CHM13#0#chrY" -> "chrY"
    """
    # Match pattern like "chr1", "chr2", "chrX", "chrY", "chr22", etc.
    # This regex matches "chr" followed by alphanumeric characters (case-insensitive)
    match = re.search(r'chr([\w]+)', path_name, re.IGNORECASE)
    if match:
        chr_full = match.group(0)  # e.g., "chr2", "chrX"
        # Normalize: lowercase 'chr', but preserve case for identifier
        chr_id = match.group(1)
        # Handle special cases like chrX, chrY, chrM - use uppercase
        if chr_id.upper() in ['X', 'Y', 'M', 'MT']:
            return f"chr{chr_id.upper()}"
        else:
            # For numeric chromosomes, use lowercase
            return f"chr{chr_id}"
    
    # Fallback: try to extract any chromosome-like identifier from parts
    parts = path_name.split('#')
    for part in parts:
        if re.match(r'^chr[\w]+$', part, re.IGNORECASE):
            # Normalize the chromosome name
            match = re.search(r'chr([\w]+)', part, re.IGNORECASE)
            if match:
                chr_id = match.group(1)
                if chr_id.upper() in ['X', 'Y', 'M', 'MT']:
                    return f"chr{chr_id.upper()}"
                else:
                    return f"chr{chr_id}"
            return part
    
    # If no chromosome found, return a sanitized version of path_name
    # Replace # with _ to make it a valid filename
    return path_name.replace('#', '_').replace('/', '_')


def get_jsonl_file_path(base_dir: Path, rank: int) -> Path:
    """
    Get the path to the JSONL file for a chunk based on its rank.
    
    Args:
        base_dir: Base directory containing rank-based subdirectories
        rank: Rank (0-based index) of the chunk in the TSV file
    
    Returns:
        Path to the JSONL file: {base_dir}/{rank}/anchors/subgraph.anchors.json.extended.jsonl
    """
    return base_dir / str(rank) / "anchors" / "subgraph.anchors.json.extended.jsonl"


def load_json(json_file: Path) -> Dict[str, Any]:
    """Load JSON file and return its contents."""
    try:
        with open(json_file, 'r') as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from {json_file}: {e}", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"Error reading {json_file}: {e}", file=sys.stderr)
        return {}


def load_jsonl(jsonl_file: Path) -> List[Dict[str, Any]]:
    """
    Load JSON/JSONL file. Handles multiple formats:
    1. JSON array file (starts with '[') - returns all items in the array
    2. JSONL file (one JSON object per line) - returns list of all objects
    3. Multi-line JSON objects - accumulates until complete
    
    Args:
        jsonl_file: Path to the JSON/JSONL file
    
    Returns:
        List of JSON objects
    """
    json_objects = []
    try:
        # First, try to read and parse as a single JSON array/object
        # This handles the case where the file is a JSON array despite .jsonl extension
        with open(jsonl_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Try parsing as a complete JSON document first
        try:
            parsed = json.loads(content)
            # If it's a list/array, return all items
            if isinstance(parsed, list):
                return parsed
            # If it's a single object, return it as a list with one item
            elif isinstance(parsed, dict):
                return [parsed]
            else:
                # Primitive value, wrap in list
                return [parsed]
        except json.JSONDecodeError:
            # Not a single JSON document, try JSONL format (one object per line)
            pass
        
        # If that failed, try JSONL format (one JSON object per line)
        # Use streaming parser to handle multi-line objects
        decoder = json.JSONDecoder()
        buffer = ""
        line_num = 0
        
        with open(jsonl_file, 'r', encoding='utf-8') as f:
            for line in f:
                line_num += 1
                buffer += line
                
                # Try to parse JSON objects from the buffer
                buffer = buffer.lstrip()
                while buffer:
                    try:
                        obj, idx = decoder.raw_decode(buffer)
                        json_objects.append(obj)
                        # Remove the parsed portion
                        buffer = buffer[idx:].lstrip()
                    except json.JSONDecodeError:
                        # Not enough data yet, continue reading
                        break
                
                # Safety check: if buffer gets too large, something might be wrong
                if len(buffer) > 100_000_000:  # 100MB limit
                    print(f"Warning: Buffer exceeded 100MB at line {line_num} of {jsonl_file}. "
                          f"Stopping to avoid memory issues.", file=sys.stderr)
                    break
        
        # Try to parse any remaining buffer
        if buffer.strip():
            try:
                obj, idx = decoder.raw_decode(buffer.strip())
                json_objects.append(obj)
            except json.JSONDecodeError:
                # If we have some objects already, warn but continue
                if json_objects:
                    print(f"Warning: Could not parse final JSON object in {jsonl_file} "
                          f"(incomplete or malformed). Processed {len(json_objects)} objects so far.", 
                          file=sys.stderr)
                else:
                    # No objects parsed at all - this is more serious
                    print(f"Error: Could not parse any JSON objects from {jsonl_file}. "
                          f"File may be malformed or in an unsupported format.", file=sys.stderr)
                                
    except FileNotFoundError:
        print(f"Error: File not found: {jsonl_file}", file=sys.stderr)
    except MemoryError:
        print(f"Error: Out of memory while reading {jsonl_file}. File may be too large.", file=sys.stderr)
    except Exception as e:
        print(f"Error reading file {jsonl_file}: {e}", file=sys.stderr)
    
    return json_objects


def combine_json_dicts(json_dicts: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Combine multiple JSON dictionaries into one.
    
    If keys overlap, the values from later dictionaries will overwrite earlier ones.
    For nested structures, we merge them appropriately.
    """
    if not json_dicts:
        return {}
    
    combined = {}
    for json_dict in json_dicts:
        if isinstance(json_dict, dict):
            combined.update(json_dict)
        else:
            # If it's not a dict, try to merge as lists or other structures
            print(f"Warning: JSON file contains non-dict structure: {type(json_dict)}", file=sys.stderr)
    
    return combined


def process_single_chromosome(chr_name: str, chr_data: Dict[str, List[Dict[str, Any]]], 
                               base_dir: Path, output_dir: Path, 
                               use_rank_based: bool, json_file_pattern: Optional[str]) -> Dict[str, Any]:
    """
    Process a single chromosome: load all JSON/JSONL files and combine them.
    
    This function is designed to be called in parallel for different chromosomes.
    
    Returns:
        Dictionary with processing results: {'chr_name': str, 'found_files': int, 
        'missing_files': List[str], 'total_json_objects': int, 'output_file': Path, 'success': bool}
    """
    result = {
        'chr_name': chr_name,
        'found_files': 0,
        'missing_files': [],
        'total_json_objects': 0,
        'output_file': None,
        'success': False
    }
    
    try:
        all_json_objects = []  # List to collect all JSON objects from all chunks
        
        # Iterate through all root chains in this chromosome
        for root_chain_id in chr_data:
            chunk_list = chr_data[root_chain_id]
            
            # Process each chunk in this root chain
            for chunk in chunk_list:
                rank = chunk['rank']
                start_node = chunk['start_node']
                end_node = chunk['end_node']
                start_orient = chunk['start_orient']
                end_orient = chunk['end_orient']
                
                if use_rank_based:
                    # Use rank-based directory structure
                    jsonl_file_path = get_jsonl_file_path(base_dir, rank)
                    
                    if not jsonl_file_path.exists():
                        result['missing_files'].append(f"rank_{rank} ({start_node}_{end_node})")
                        continue
                    
                    # Load JSONL file (contains multiple JSON objects, one per line)
                    json_objects = load_jsonl(jsonl_file_path)
                    if json_objects:
                        all_json_objects.extend(json_objects)
                        result['total_json_objects'] += len(json_objects)
                        result['found_files'] += 1
                else:
                    # Use pattern-based file finding
                    if not json_file_pattern:
                        result['missing_files'].append(f"{start_node}_{end_node}")
                        continue
                    
                    json_file_path = base_dir / json_file_pattern.format(
                        start_node=start_node,
                        end_node=end_node,
                        start_orient=start_orient,
                        end_orient=end_orient
                    )
                    if not json_file_path.exists():
                        result['missing_files'].append(f"{start_node}_{end_node}")
                        continue
                    
                    # Load JSON file
                    json_data = load_json(json_file_path)
                    if json_data:
                        all_json_objects.append(json_data)
                        result['found_files'] += 1
                        result['total_json_objects'] += 1
        
        if not all_json_objects:
            print(f"  Error: No JSON objects found for {chr_name}. Skipping.", file=sys.stderr)
            return result
        
        # Combine JSON objects
        # For JSONL files (rank-based), preserve all objects as a list since each line is a separate record
        # For regular JSON files, combine dicts if they're all dicts
        if use_rank_based:
            # JSONL format: output as a list of all JSON objects
            combined_json = all_json_objects
        elif all_json_objects and all(isinstance(obj, dict) for obj in all_json_objects):
            # Regular JSON files: combine dicts (later ones overwrite earlier ones for duplicate keys)
            combined_json = combine_json_dicts(all_json_objects)
        else:
            # Mixed types: store as a list wrapped in a dict with metadata
            combined_json = {'objects': all_json_objects, 'count': len(all_json_objects)}
        
        # Write combined JSON file
        output_file = output_dir / f"{chr_name}_combined.json"
        with open(output_file, 'w') as f:
            json.dump(combined_json, f, indent=2)
        
        result['output_file'] = output_file
        result['success'] = True
        
    except Exception as e:
        print(f"  Error processing chromosome {chr_name}: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
    
    return result


def process_chunks(tsv_file: Path, base_dir: Path, output_dir: Path, 
                   use_rank_based: bool = True, json_file_pattern: Optional[str] = None,
                   num_workers: Optional[int] = None) -> None:
    """
    Process chunks from TSV file and combine JSON/JSONL files per chromosome.
    
    Args:
        tsv_file: Path to TSV file with chunk information
        base_dir: Base directory containing rank-based subdirectories (if use_rank_based=True)
                  or directory containing JSON files (if use_rank_based=False)
        output_dir: Directory to write combined JSON files
        use_rank_based: If True, use rank-based directory structure ({base_dir}/{rank}/anchors/...)
                        If False, use traditional file finding by node IDs
        json_file_pattern: Optional pattern for JSON file naming (e.g., "{start_node}_{end_node}.json")
                          Only used if use_rank_based=False
        num_workers: Number of parallel workers to use. If None, uses CPU count.
    """
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read TSV file and group by chromosome
    # Store rank (0-based index from TSV line number) with each chunk
    chr_chunks = {}  # {chr1: {c1234567890-1234567891: [{rank: 0, start_node: ..., end_node: ..., ...}, ...]}}
    
    with open(tsv_file, 'r') as f:
        header = f.readline().strip().split('\t')
        
        # Find column indices
        try:
            start_node_idx = header.index('start_node')
            start_orient_idx = header.index('start_orientation')
            end_node_idx = header.index('end_node')
            end_orient_idx = header.index('end_orientation')
            path_name_idx = header.index('path_name')
            root_chain_id_idx = header.index('root_chain_id')
        
        except ValueError as e:
            print(f"Error: Required column not found in TSV file: {e}", file=sys.stderr)
            print(f"Available columns: {header}", file=sys.stderr)
            sys.exit(1)
        
        # Track rank (0-based index: first data row is rank 0)
        rank = 0
        for line_num, line in enumerate(f, start=2):
            if not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < len(header):
                print(f"Warning: Line {line_num} has fewer fields than header. Skipping.", file=sys.stderr)
                continue
            
            try:
                start_node = int(fields[start_node_idx])
                start_orient = fields[start_orient_idx]
                end_node = int(fields[end_node_idx])
                end_orient = fields[end_orient_idx]
                path_name = fields[path_name_idx]
                root_chain_id = fields[root_chain_id_idx].strip()
                
                # Use a default root_chain_id if empty
                if not root_chain_id:
                    root_chain_id = f"c{start_node}-{end_node}"

                # Extract chromosome name
                chr_name = extract_chr_name(path_name)
                
                # Store chunk information with rank
                if chr_name not in chr_chunks:
                    chr_chunks[chr_name] = {}
                if root_chain_id not in chr_chunks[chr_name]:
                    chr_chunks[chr_name][root_chain_id] = []
                chr_chunks[chr_name][root_chain_id].append({
                    'rank': rank,
                    'start_node': start_node,
                    'end_node': end_node,
                    'start_orient': start_orient,
                    'end_orient': end_orient
                })
                
                # Increment rank for next chunk
                rank += 1
            except (ValueError, IndexError) as e:
                print(f"Warning: Error parsing line {line_num}: {e}. Skipping.", file=sys.stderr)
                continue
    
    print(f"Found {len(chr_chunks)} chromosomes", file=sys.stderr)
    for chr_name in list(chr_chunks.keys()):
        root_chains_in_chr = list(chr_chunks[chr_name].keys())
        print(f"  {chr_name}: {len(root_chains_in_chr)} root chains", file=sys.stderr)
    
    # Determine number of workers
    if num_workers is None:
        num_workers = multiprocessing.cpu_count()
    num_workers = min(num_workers, len(chr_chunks))  # Don't use more workers than chromosomes
    
    print(f"\nProcessing {len(chr_chunks)} chromosomes using {num_workers} workers...", file=sys.stderr)
    
    # Process chromosomes in parallel
    if num_workers > 1 and len(chr_chunks) > 1:
        # Use parallel processing
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            # Submit all chromosome processing tasks
            future_to_chr = {
                executor.submit(
                    process_single_chromosome,
                    chr_name,
                    chr_chunks[chr_name],
                    base_dir,
                    output_dir,
                    use_rank_based,
                    json_file_pattern
                ): chr_name
                for chr_name in chr_chunks
            }
            
            # Collect results as they complete
            completed = 0
            for future in as_completed(future_to_chr):
                chr_name = future_to_chr[future]
                completed += 1
                try:
                    result = future.result()
                    if result['success']:
                        print(f"[{completed}/{len(chr_chunks)}] Completed {chr_name}: "
                              f"{result['found_files']} files, {result['total_json_objects']} JSON objects -> "
                              f"{result['output_file'].name}", file=sys.stderr)
                        if result['missing_files']:
                            print(f"  Warning: {len(result['missing_files'])} missing files: "
                                  f"{result['missing_files'][:5]}...", file=sys.stderr)
                    else:
                        print(f"[{completed}/{len(chr_chunks)}] Failed {chr_name}: "
                              f"No JSON objects found", file=sys.stderr)
                except Exception as e:
                    print(f"[{completed}/{len(chr_chunks)}] Error processing {chr_name}: {e}", file=sys.stderr)
    else:
        # Sequential processing (single worker or single chromosome)
        for chr_name in chr_chunks:
            print(f"\nProcessing chromosome {chr_name}...", file=sys.stderr)
            result = process_single_chromosome(
                chr_name,
                chr_chunks[chr_name],
                base_dir,
                output_dir,
                use_rank_based,
                json_file_pattern
            )
            if result['success']:
                print(f"  Combined {result['found_files']} files ({result['total_json_objects']} JSON objects) "
                      f"into {result['output_file']}", file=sys.stderr)
                if isinstance(result.get('output_file'), Path):
                    # Load to check structure for reporting
                    try:
                        with open(result['output_file'], 'r') as f:
                            combined_json = json.load(f)
                        if isinstance(combined_json, list):
                            print(f"  Total JSON objects in output: {len(combined_json)}", file=sys.stderr)
                        elif isinstance(combined_json, dict) and 'objects' not in combined_json:
                            print(f"  Total keys in combined JSON: {len(combined_json)}", file=sys.stderr)
                    except:
                        pass
                if result['missing_files']:
                    print(f"  Warning: Could not find JSON/JSONL files for {len(result['missing_files'])} chunks: "
                          f"{result['missing_files'][:10]}...", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Combine JSON/JSONL files per chromosome based on chunk information from TSV file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use rank-based directory structure (default) with parallel processing
  %(prog)s -t chunks.tsv -b /path/to/base/dir -o output/ -j 8

  # Use traditional file finding by node IDs
  %(prog)s -t chunks.tsv -b /path/to/json/dir -o output/ --no-rank-based

  # Use custom file pattern with 4 workers
  %(prog)s -t chunks.tsv -b /path/to/json/dir -o output/ --no-rank-based -p "{start_node}_{end_node}.json" -j 4
        """
    )
    parser.add_argument('-t', '--tsv-file', type=Path, required=True,
                       help='Path to TSV file with chunk information (must have start_node, end_node, path_name, root_chain_id columns)')
    parser.add_argument('-b', '--base-dir', type=Path, required=True,
                       help='Base directory containing rank-based subdirectories (if --use-rank-based) or JSON files (if --no-rank-based)')
    parser.add_argument('-o', '--output-dir', type=Path, default=Path('combined_jsons'),
                       help='Directory to write combined JSON files (default: combined_jsons)')
    parser.add_argument('--use-rank-based', action='store_true', default=True,
                       help='Use rank-based directory structure: {base_dir}/{rank}/anchors/subgraph.anchors.json.extended.jsonl (default: True)')
    parser.add_argument('--no-rank-based', dest='use_rank_based', action='store_false',
                       help='Use traditional file finding by node IDs instead of rank-based structure')
    parser.add_argument('-p', '--pattern', default=None,
                       help='Pattern for JSON file naming (e.g., "{start_node}_{end_node}.json"). Required when using --no-rank-based')
    parser.add_argument('-j', '--jobs', type=int, default=None,
                       help='Number of parallel workers to use (default: number of CPU cores)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.tsv_file.exists():
        print(f"Error: TSV file not found: {args.tsv_file}", file=sys.stderr)
        sys.exit(1)
    
    if not args.base_dir.exists():
        print(f"Error: Base directory not found: {args.base_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Process chunks
    process_chunks(args.tsv_file, args.base_dir, args.output_dir, 
                  use_rank_based=args.use_rank_based, json_file_pattern=args.pattern,
                  num_workers=args.jobs)


if __name__ == '__main__':
    main()

