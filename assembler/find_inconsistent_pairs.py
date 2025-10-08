import argparse
import json
import os
from collections import defaultdict, Counter

def parse_anchors_json(file_path):
    """
    Parses the anchors.json file and builds the anchor-to-reads mapping.
    The format is {anchor_id: {read_id: strand}}.
    """
    anchor_to_reads = defaultdict(dict)
    with open(file_path, 'r') as f:
        data = json.load(f)
        
        for anchor_data in data:
            if not isinstance(anchor_data, list) or len(anchor_data) != 2:
                continue
            
            anchor_id = anchor_data[0]
            reads_list = anchor_data[1]
            
            for read_info in reads_list:
                if not isinstance(read_info, list) or len(read_info) < 2:
                    continue
                read_id = read_info[0]
                strand = read_info[1]
                anchor_to_reads[anchor_id][read_id] = strand
    return anchor_to_reads

def build_orientation_dict(anchor_to_reads):
    """
    Builds a dictionary of relative orientations for each read pair on shared anchors.
    
    Returns:
        Dict[Tuple[read_A, read_B], Dict[anchor_id, orientation]]
    """
    read_pair_orientations = defaultdict(dict)
    
    for anchor_id, reads in anchor_to_reads.items():
        if len(reads) < 2:
            continue
        
        read_ids = sorted(list(reads.keys()))
        for i in range(len(read_ids)):
            for j in range(i + 1, len(read_ids)):
                read_A, read_B = read_ids[i], read_ids[j]
                strand_A = reads[read_A]
                strand_B = reads[read_B]
                
                orientation = "same" if strand_A == strand_B else "different"
                read_pair = (read_A, read_B)
                read_pair_orientations[read_pair][anchor_id] = orientation
                
    return read_pair_orientations

def find_conflicting_pairs(read_pair_orientations):
    """
    Identifies read pairs that have conflicting orientations across different anchors.
    
    Returns:
        A dictionary of conflicting pairs and their orientation data.
    """
    conflicts = {}
    for read_pair, anchor_orientations in read_pair_orientations.items():
        orientations = anchor_orientations.values()
        if "same" in orientations and "different" in orientations:
            conflicts[read_pair] = anchor_orientations
    return conflicts

def tally_inconsistent_reads(conflicting_pairs):
    """
    Counts the occurrences of each read in the set of conflicting pairs.
    
    Returns:
        A list of (read_id, count) tuples, sorted by count descending.
    """
    read_counts = Counter()
    for read_A, read_B in conflicting_pairs.keys():
        read_counts[read_A] += 1
        read_counts[read_B] += 1
    return sorted(read_counts.items(), key=lambda item: item[1], reverse=True)

def write_pruned_anchors(input_file, output_dir, reads_to_remove):
    """
    Reads an anchors.json file, removes specified reads, and writes a new pruned version.
    """
    with open(input_file, 'r') as f:
        original_data = json.load(f)
    
    pruned_data = []
    for anchor_data in original_data:
        if not isinstance(anchor_data, list) or len(anchor_data) != 2:
            pruned_data.append(anchor_data)
            continue
            
        anchor_id, reads_list = anchor_data
        
        # Filter out reads that are in the removal set
        pruned_reads_list = [
            read_info for read_info in reads_list if read_info[0] not in reads_to_remove
        ]
        
        pruned_data.append([anchor_id, pruned_reads_list])
        
    output_path = os.path.join(output_dir, "anchors.pruned.json")
    with open(output_path, 'w') as f:
        json.dump(pruned_data, f, indent=4)
        
    print(f"Pruned anchors file saved to: {output_path}")

def main():
    """Main function to run the analysis script."""
    parser = argparse.ArgumentParser(
        description="Finds inconsistent read pairs from an anchors.json file and reports problematic reads."
    )
    parser.add_argument(
        "anchors_file", 
        help="Path to the anchors.json file."
    )
    parser.add_argument(
        "output_dir",
        help="Path to the directory where output files will be saved."
    )
    args = parser.parse_args()

    print(f"Analyzing file: {args.anchors_file}")
    
    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    anchor_map = parse_anchors_json(args.anchors_file)
    orientation_dict = build_orientation_dict(anchor_map)
    conflicting_pairs = find_conflicting_pairs(orientation_dict)
    
    # --- Write Inconsistent Pairs Report ---
    pairs_report_path = os.path.join(args.output_dir, "inconsistent_pairs.txt")
    with open(pairs_report_path, 'w') as f:
        if conflicting_pairs:
            num_conflicts = len(conflicting_pairs)
            f.write(f"Found {num_conflicts} read pairs with conflicting orientations.\n")
            f.write("\n--- Inconsistent Read Pairs Report ---\n")
            
            sorted_pairs = sorted(conflicting_pairs.items(), key=lambda item: item[0][0])
            
            for pair, orientations in sorted_pairs:
                f.write(f"\n- Pair {pair} has conflicting orientations:\n")
                sorted_anchors = sorted(orientations.items(), key=lambda item: item[0])
                for anchor_id, orientation in sorted_anchors:
                    f.write(f"  - Anchor '{anchor_id}': {orientation}\n")
            f.write("\n--------------------------------------\n")
        else:
            f.write("No read pairs with conflicting orientations were found.\n")
    
    print(f"Inconsistent pairs report saved to: {pairs_report_path}")

    # --- Write Problematic Reads Report ---
    if conflicting_pairs:
        reads_report_path = os.path.join(args.output_dir, "problematic_reads_report.txt")
        inconsistent_read_counts = tally_inconsistent_reads(conflicting_pairs)
        
        with open(reads_report_path, 'w') as f:
            f.write("Read_ID\tFrequency_in_Conflicts\n")
            for read_id, count in inconsistent_read_counts:
                f.write(f"{read_id}\t{count}\n")
        
        print(f"Problematic reads report saved to: {reads_report_path}")

        # --- Prune and Write New anchors.json ---
        reads_to_remove = {
            read_id for read_id, count in inconsistent_read_counts if count > 1
        }
        
        if reads_to_remove:
            print(f"\nFound {len(reads_to_remove)} reads with frequency > 1. Pruning them from anchors.json.")
            write_pruned_anchors(args.anchors_file, args.output_dir, reads_to_remove)
        else:
            print("\nNo reads found with conflict frequency > 1. No pruned anchors file will be written.")

if __name__ == "__main__":
    main()
