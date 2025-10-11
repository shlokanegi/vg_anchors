import argparse
import csv
import sys
import json

def find_mixed_strand_rows(file_path):
    """
    Parses a journeys.csv file to find reads that contain both '+' and '-' anchors.
    
    Args:
        file_path (str): The path to the journeys.csv file.
        
    Returns:
        A list of rows (as lists of strings) for reads that have mixed anchor strands.
    """
    mixed_strand_rows = []
    with open(file_path, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            
            anchors_journey = row[1:]
            
            has_plus = False
            has_minus = False
            
            for anchor_entry in anchors_journey:
                if not anchor_entry:
                    continue
                
                if anchor_entry.endswith('+'):
                    has_plus = True
                elif anchor_entry.endswith('-'):
                    has_minus = True
                
                # Optimization: if we've found both, we can stop checking this read
                if has_plus and has_minus:
                    mixed_strand_rows.append(row)
                    break  # Move to the next row
                    
    return mixed_strand_rows

def load_read_mapping(file_path):
    """
    Loads a read ID to read name mapping from a CSV file.
    The CSV is expected to have a header. The mapping is from the second column
    (OrientedReadId) to the third column (ReadName).
    
    Args:
        file_path (str): Path to the CSV file.
        
    Returns:
        A dictionary mapping read IDs to read names.
    """
    read_map = {}
    with open(file_path, 'r') as f:
        reader = csv.reader(f)
        try:
            next(reader)  # Skip header
        except StopIteration:
            return {}  # Handle empty file
        for row in reader:
            if len(row) >= 3:
                read_map[row[1]] = row[2]
    return read_map

def prune_anchors(anchors_data, bad_read_names):
    """
    Removes reads with mixed strands from the anchors data.
    
    Args:
        anchors_data (list): The original anchors data from anchors.json.
        bad_read_names (set): A set of read names to be removed.
        
    Returns:
        A new list of anchors with the bad reads removed.
    """
    pruned_anchors = []
    for anchor_info in anchors_data:
        anchor_id = anchor_info[0]
        reads = anchor_info[1]
        
        pruned_reads = [read for read in reads if read[0] not in bad_read_names]
        
        if pruned_reads:
            pruned_anchors.append([anchor_id, pruned_reads])
            
    return pruned_anchors

def main():
    """Main function to run the analysis script."""
    parser = argparse.ArgumentParser(
        description="Outputs journey data for reads from a Journeys.csv file that contain anchors with both '+' and '-' orientations. The output is in CSV format."
    )
    parser.add_argument(
        "journeys_file", 
        help="Path to the Journeys.csv file."
    )
    parser.add_argument(
        "-o", "--output-prefix",
        help="Path to the output prefix.",
        default="mixed_strand_reads"
    )
    parser.add_argument(
        "-a", "--anchors-json",
        help="Path to the anchors.json file. If provided, the script will generate a pruned version of this file."
    )
    parser.add_argument(
        "-r", "--read-mapping-file",
        help="Path to a CSV file mapping read IDs to read names. Required if --anchors-json is provided."
    )
    args = parser.parse_args()

    if args.anchors_json and not args.read_mapping_file:
        parser.error("--read-mapping-file is required when --anchors-json is provided.")

    # Print status messages to stderr to keep stdout clean for redirection
    print(f"Analyzing file: {args.journeys_file}", file=sys.stderr)
    
    # Run the analysis to find rows for reads with mixed strands
    bad_read_rows = find_mixed_strand_rows(args.journeys_file)
    
    if bad_read_rows:
        # Sort rows numerically based on the read_id (e.g., '10-1' > '2-1')
        try:
            sorted_rows = sorted(
                bad_read_rows, 
                key=lambda r: (int(r[0].split('-')[0]), int(r[0].split('-')[1]))
            )
        except (ValueError, IndexError):
            # Fallback to simple string sort if read_id format is unexpected
            sorted_rows = sorted(bad_read_rows, key=lambda r: r[0])

        # Output another file with just the switch anchors, i.e. where the read switches strands
        # Emit pairs as: [read_id, anchor_i, anchor_{i+1}]
        switch_anchors_rows = []
        for row in sorted_rows:
            read_id = row[0]
            anchors_journey = row[1:]
            for i in range(len(anchors_journey) - 1):
                a, b = anchors_journey[i], anchors_journey[i + 1]
                if not a or not b:
                    continue
                if (a.endswith('+') and b.endswith('-')) or (a.endswith('-') and b.endswith('+')):
                    switch_anchors_rows.append([read_id, a, b])

        # Write outputs with proper newline handling
        with open(args.output_prefix + ".mixed_strand_reads.csv", "w", newline="") as f1:
            writer = csv.writer(f1)
            writer.writerows(sorted_rows)

        with open(args.output_prefix + ".switch_anchors.csv", "w", newline="") as f2:
            writer2 = csv.writer(f2)
            writer2.writerows(switch_anchors_rows)

        if args.anchors_json:
            print(f"Pruning anchors from {args.anchors_json}", file=sys.stderr)
            
            bad_read_ids = {row[0] for row in bad_read_rows}
            
            read_id_to_name = load_read_mapping(args.read_mapping_file)
            
            bad_read_names = {read_id_to_name[read_id] for read_id in bad_read_ids if read_id in read_id_to_name}
            
            print(f"Found {len(bad_read_ids)} mixed-strand read IDs, mapping to {len(bad_read_names)} unique read names to be pruned.", file=sys.stderr)

            with open(args.anchors_json, 'r') as f:
                anchors_data = json.load(f)
                
            pruned_anchors_data = prune_anchors(anchors_data, bad_read_names)
            
            if args.anchors_json.endswith('.json'):
                output_json_path = args.anchors_json[:-5] + ".strand-pruned.json"
            else:
                output_json_path = args.anchors_json + ".strand-pruned.json"
                
            with open(output_json_path, 'w') as f:
                json.dump(pruned_anchors_data, f, indent=2)
                
            print(f"Wrote pruned anchors to {output_json_path}", file=sys.stderr)
    
    else:
        print("No reads with mixed anchor strands were found.", file=sys.stderr)

if __name__ == "__main__":
    main()
