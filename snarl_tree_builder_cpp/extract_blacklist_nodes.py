#!/usr/bin/env python3
"""
Step 4: Extract Blacklist Nodes from Centromere Regions

This script:
1. Reads censat_merged_regions.bed
2. For each region, runs gbz-query to extract the subgraph
3. Parses the GFA output to extract node IDs from 'S' lines
4. Combines nodes per chromosome (handling multiple regions per chromosome)
5. Saves unique, sorted node IDs as binary int64_t files (one per chromosome)

Usage:
    python3 extract_blacklist_nodes.py [--threads N] [--gbz-db PATH] [--sample SAMPLE] [--input-bed PATH] [--output-dir PATH]

Options:
    --threads N: Number of threads to use (default: 16)
    --gbz-db PATH: Path to GBZ database (default: /data/tmp/PAW70337-16-sampled.gbz.db)
    --sample SAMPLE: Sample name (default: CHM13)
    --input-bed PATH: Path to input BED file (default: /data/tmp/censat_merged_regions.bed)
    --output-dir PATH: Path to output directory (default: /data/tmp/blacklist_nodes)
"""

import os
import sys
import argparse
import subprocess
import struct
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
import time

# Configuration
DEFAULT_GBZ_DB = "/data/tmp/PAW70337-16-sampled.gbz.db"
# DEFAULT_GBZ_DB = "/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/results_hs/graph/PAW70337/gbz_db/PAW70337-16-sampled.gbz.db"
DEFAULT_SAMPLE = "CHM13"
DEFAULT_THREADS = 16
INPUT_BED = "/data/tmp/censat_merged_regions.bed"
OUTPUT_DIR = "/data/tmp/blacklist_nodes"

# Chromosomes to skip (chrY causes issues)
SKIP_CHROMOSOMES = {"chrY"}

# Lock for thread-safe operations
print_lock = Lock()
nodes_lock = Lock()


def log(message):
    """Thread-safe logging."""
    with print_lock:
        print(message, flush=True)


def parse_gfa_nodes(gfa_content):
    """
    Parse GFA content and extract node IDs from 'S' (segment) lines.
    
    GFA format for S lines: S<tab>node_id<tab>sequence[<tab>optional_tags]
    
    Returns a set of node IDs (as int64).
    """
    nodes = set()
    for line in gfa_content.split('\n'):
        if line.startswith('S\t'):
            parts = line.split('\t')
            if len(parts) >= 2:
                try:
                    node_id = int(parts[1])
                    nodes.add(node_id)
                except ValueError:
                    # Skip if node_id is not a valid integer
                    pass
    return nodes


def run_gbz_query(chrom, start, end, gbz_db, sample):
    """
    Run gbz-query for a specific region and return extracted node IDs.
    
    Args:
        chrom: Chromosome name (e.g., 'chr1')
        start: Start position (0-based)
        end: End position
        gbz_db: Path to GBZ database
        sample: Sample name (e.g., 'CHM13')
    
    Returns:
        tuple: (chrom, set of node IDs, error message or None)
    """
    interval = f"{start}..{end}"
    
    cmd = [
        "query",
        "--sample", sample,
        "--contig", chrom,
        "--interval", interval,
        "--snarls", gbz_db
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=1800  # 30 minute timeout per region
        )
        
        if result.returncode != 0:
            error_msg = result.stderr.strip() if result.stderr else f"Exit code {result.returncode}"
            return (chrom, set(), f"gbz-query failed: {error_msg}")
        
        # Parse GFA from stdout
        nodes = parse_gfa_nodes(result.stdout)
        return (chrom, nodes, None)
        
    except subprocess.TimeoutExpired:
        return (chrom, set(), "Timeout (>10 min)")
    except Exception as e:
        return (chrom, set(), str(e))


def process_region(region_info, gbz_db, sample):
    """
    Process a single centromere region.
    
    Args:
        region_info: tuple of (chrom, start, end, region_idx, total_regions)
        gbz_db: Path to GBZ database
        sample: Sample name
    
    Returns:
        tuple: (chrom, set of node IDs, error message or None, region_info)
    """
    chrom, start, end, region_idx, total_regions = region_info
    
    log(f"  [{region_idx}/{total_regions}] Processing {chrom}:{start}-{end} ({end-start:,} bp)...")
    
    start_time = time.time()
    chrom_result, nodes, error = run_gbz_query(chrom, start, end, gbz_db, sample)
    elapsed = time.time() - start_time
    
    if error:
        log(f"  [{region_idx}/{total_regions}] ERROR {chrom}:{start}-{end}: {error}")
    else:
        log(f"  [{region_idx}/{total_regions}] Done {chrom}:{start}-{end}: {len(nodes):,} nodes in {elapsed:.1f}s")
    
    return (chrom, nodes, error, (chrom, start, end))


def read_bed_file(bed_path):
    """
    Read the BED file and return list of regions.
    
    Returns:
        list of tuples: (chrom, start, end)
    """
    regions = []
    
    with open(bed_path, 'r') as f:
        for line in f:
            # Skip header/comment lines
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            
            # Skip excluded chromosomes
            if chrom in SKIP_CHROMOSOMES:
                continue
            
            regions.append((chrom, start, end))
    
    return regions


def save_nodes_binary(nodes, output_path):
    """
    Save sorted node IDs as binary int64_t file.
    
    Args:
        nodes: set or list of node IDs
        output_path: Path to output file
    """
    sorted_nodes = sorted(nodes)
    
    with open(output_path, 'wb') as f:
        for node_id in sorted_nodes:
            # Pack as signed 64-bit integer (little-endian)
            f.write(struct.pack('<q', node_id))
    
    return len(sorted_nodes)


def main():
    parser = argparse.ArgumentParser(
        description="Extract blacklist nodes from centromere regions using gbz-query"
    )
    parser.add_argument(
        "--threads", "-t",
        type=int,
        default=DEFAULT_THREADS,
        help=f"Number of threads for parallel processing (default: {DEFAULT_THREADS})"
    )
    parser.add_argument(
        "--gbz-db",
        type=str,
        default=DEFAULT_GBZ_DB,
        help=f"Path to GBZ database (default: {DEFAULT_GBZ_DB})"
    )
    parser.add_argument(
        "--sample",
        type=str,
        default=DEFAULT_SAMPLE,
        help=f"Sample name for gbz-query (default: {DEFAULT_SAMPLE})"
    )
    parser.add_argument(
        "--input-bed",
        type=str,
        default=INPUT_BED,
        help=f"Input BED file with centromere regions (default: {INPUT_BED})"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=OUTPUT_DIR,
        help=f"Output directory for binary node files (default: {OUTPUT_DIR})"
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.input_bed):
        print(f"ERROR: Input BED file not found: {args.input_bed}")
        sys.exit(1)
    
    if not os.path.exists(args.gbz_db):
        print(f"ERROR: GBZ database not found: {args.gbz_db}")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("=" * 70)
    print("Step 4: Extract Blacklist Nodes from Centromere Regions")
    print("=" * 70)
    print(f"Input BED:    {args.input_bed}")
    print(f"GBZ database: {args.gbz_db}")
    print(f"Sample:       {args.sample}")
    print(f"Threads:      {args.threads}")
    print(f"Output dir:   {args.output_dir}")
    print(f"Skipping:     {', '.join(SKIP_CHROMOSOMES)}")
    print("=" * 70)
    
    # Read regions from BED file
    print("\n[1] Reading centromere regions from BED file...")
    regions = read_bed_file(args.input_bed)
    print(f"    Found {len(regions)} regions (after filtering)")
    
    # Count regions per chromosome
    chrom_counts = defaultdict(int)
    for chrom, _, _ in regions:
        chrom_counts[chrom] += 1
    print(f"    Chromosomes: {len(chrom_counts)}")
    for chrom in sorted(chrom_counts.keys()):
        print(f"      {chrom}: {chrom_counts[chrom]} regions")
    
    # Prepare region info with indices
    region_infos = [
        (chrom, start, end, i+1, len(regions))
        for i, (chrom, start, end) in enumerate(regions)
    ]
    
    # Process regions in parallel
    print(f"\n[2] Extracting subgraphs using gbz-query ({args.threads} threads)...")
    
    # Collect nodes per chromosome
    chromosome_nodes = defaultdict(set)
    failed_regions = []
    
    start_time = time.time()
    
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        # Submit all tasks
        futures = {
            executor.submit(process_region, info, args.gbz_db, args.sample): info
            for info in region_infos
        }
        
        # Collect results as they complete
        for future in as_completed(futures):
            chrom, nodes, error, region_info = future.result()
            
            if error:
                failed_regions.append((region_info, error))
            else:
                with nodes_lock:
                    chromosome_nodes[chrom].update(nodes)
    
    elapsed = time.time() - start_time
    print(f"\n    Total extraction time: {elapsed:.1f}s")
    
    # Report failures
    if failed_regions:
        print(f"\n    WARNING: {len(failed_regions)} regions failed:")
        for (chrom, start, end), error in failed_regions:
            print(f"      {chrom}:{start}-{end}: {error}")
    
    # Save nodes to binary files
    print(f"\n[3] Saving nodes to binary files...")
    
    total_nodes = 0
    for chrom in sorted(chromosome_nodes.keys()):
        nodes = chromosome_nodes[chrom]
        output_path = os.path.join(args.output_dir, f"{chrom}_blacklist_nodes.bin")
        
        num_nodes = save_nodes_binary(nodes, output_path)
        total_nodes += num_nodes
        
        print(f"    {chrom}: {num_nodes:,} unique nodes -> {output_path}")
    
    # Also save a combined file with all nodes
    all_nodes = set()
    for nodes in chromosome_nodes.values():
        all_nodes.update(nodes)
    
    combined_path = os.path.join(args.output_dir, "all_blacklist_nodes.bin")
    combined_count = save_nodes_binary(all_nodes, combined_path)
    print(f"\n    Combined: {combined_count:,} unique nodes -> {combined_path}")
    
    # Write a summary TSV
    summary_path = os.path.join(args.output_dir, "blacklist_summary.tsv")
    with open(summary_path, 'w') as f:
        f.write("chromosome\tnum_regions\tnum_nodes\tbinary_file\n")
        for chrom in sorted(chromosome_nodes.keys()):
            f.write(f"{chrom}\t{chrom_counts[chrom]}\t{len(chromosome_nodes[chrom])}\t{chrom}_blacklist_nodes.bin\n")
        f.write(f"ALL\t{len(regions)}\t{combined_count}\tall_blacklist_nodes.bin\n")
    print(f"    Summary: {summary_path}")
    
    # Write failed regions if any
    if failed_regions:
        failed_path = os.path.join(args.output_dir, "failed_regions.tsv")
        with open(failed_path, 'w') as f:
            f.write("chromosome\tstart\tend\terror\n")
            for (chrom, start, end), error in failed_regions:
                f.write(f"{chrom}\t{start}\t{end}\t{error}\n")
        print(f"    Failed regions: {failed_path}")
    
    print("\n" + "=" * 70)
    print("COMPLETE")
    print(f"  Regions processed: {len(regions)}")
    print(f"  Regions failed:    {len(failed_regions)}")
    print(f"  Chromosomes:       {len(chromosome_nodes)}")
    print(f"  Total unique nodes: {combined_count:,}")
    print("=" * 70)

    # Move OUTPUT_DIR to /private/groups/migalab/shnegi/vg_anchors_project/vg_anchors/snarl_tree_builder_cpp/blacklist_nodes
    shutil.move(args.output_dir, "/private/groups/migalab/shnegi/vg_anchors_project/vg_anchors/snarl_tree_builder_cpp/blacklist_nodes")

if __name__ == "__main__":
    main()


