#!/usr/bin/env python3
# This script prepares data for visualization by transforming raw anchor data,
# mapping anchors to snarls, and calculating anchor connectivity.

import json
import csv
import sys
import argparse
import itertools
from collections import defaultdict

def transform_raw_json(raw_input_data):
    """
    Transforms the initial list-based JSON into a dictionary format.
    Input format: [["anchor_id1", [[read1], [read2]]], ["anchor_id2", [[read3]]]]
    Output format: {"anchor_id1": [{"read_id":...}, ...], "anchor_id2": [...]}
    """
    output_dict = {}
    for item in raw_input_data:
        if not (isinstance(item, list) and len(item) == 2):
            print(f"Warning: Skipping malformed top-level item in input JSON: {item}", file=sys.stderr)
            continue

        anchor_id, reads = item
        if anchor_id not in output_dict:
            output_dict[anchor_id] = []
        
        for read in reads:
            if not (isinstance(read, list) and len(read) == 4):
                print(f"Warning: Skipping malformed read data for anchor {anchor_id}: {read}", file=sys.stderr)
                continue
            
            read_info_dict = {
                "read_id": read[0],
                "read_strand": read[1],
                "read_start": read[2],
                "read_end": read[3]
            }
            output_dict[anchor_id].append(read_info_dict)
            
    return output_dict

def create_snarl_to_anchor_mapping(transformed_data, tsv_filename):
    """
    Uses a TSV file to map anchor IDs from the transformed data to snarl IDs.
    """
    relevant_anchor_ids = set(transformed_data.keys())
    snarl_to_anchors_map = defaultdict(list)

    try:
        with open(tsv_filename, 'r', newline='') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)  # Read and store the header row
            
            anchor_col_idx = header.index('Anchor_path')
            snarl_col_idx = header.index('snarl_id')

            for row in reader:
                if len(row) > max(anchor_col_idx, snarl_col_idx):
                    anchor_id = row[anchor_col_idx]
                    snarl_id = row[snarl_col_idx]
                    if anchor_id in relevant_anchor_ids and anchor_id not in snarl_to_anchors_map[snarl_id]:
                        snarl_to_anchors_map[snarl_id].append(anchor_id)
    except (FileNotFoundError, ValueError, IndexError) as e:
        print(f"Error processing TSV file '{tsv_filename}': {e}", file=sys.stderr)
        sys.exit(1)

    return dict(snarl_to_anchors_map)

def create_anchor_connectivity(transformed_data):
    """
    Calculates how many reads are shared between each pair of anchors.
    """
    read_to_anchors = defaultdict(set)
    for anchor_id, reads in transformed_data.items():
        for read in reads:
            read_to_anchors[read['read_id']].add(anchor_id)

    connectivity = defaultdict(lambda: defaultdict(int))
    for anchors_set in read_to_anchors.values():
        if len(anchors_set) > 1:
            for anchor1, anchor2 in itertools.combinations(sorted(list(anchors_set)), 2):
                connectivity[anchor1][anchor2] += 1
    
    final_connectivity = defaultdict(dict)
    for anchor1, connections in connectivity.items():
        for anchor2, count in connections.items():
            final_connectivity[anchor1][anchor2] = {"common_read_count": count}

    return dict(final_connectivity)

def create_read_metadata(transformed_data):
    """
    Generates metadata for each unique read, including its journey across anchors,
    ordered by the read's start position in each anchor.
    """
    # Use a temporary structure to collect all anchor visits for each read
    read_journeys = defaultdict(list)
    read_strands = {}

    for anchor_id, reads in transformed_data.items():
        for read in reads:
            read_id = read['read_id']
            
            # Store the strand once, assuming it's consistent for a given read
            if read_id not in read_strands:
                read_strands[read_id] = "+" if read['read_strand'] == '0' else "-"
            
            # Store the anchor visit as a tuple of (start_position, anchor_id)
            read_journeys[read_id].append((read['read_start'], anchor_id))

    # Process the collected data to build the final metadata dictionary
    final_metadata = {}
    for read_id, journey_tuples in read_journeys.items():
        # Sort the journey by the start position (the first element of the tuple)
        sorted_journey = sorted(journey_tuples)
        
        # Extract the anchor IDs from the sorted list, ensuring uniqueness while preserving order
        final_journey_path = []
        seen_anchors = set()
        for _, anchor_id in sorted_journey:
            if anchor_id not in seen_anchors:
                final_journey_path.append(anchor_id)
                seen_anchors.add(anchor_id)
        
        final_metadata[read_id] = {
            "read_strand": read_strands[read_id],
            "journey": final_journey_path
        }
        
    return final_metadata

def main():
    """
    Main function to parse arguments and run the data preparation workflow.
    """
    parser = argparse.ArgumentParser(description="Transform raw anchor data and generate snarl, connectivity, and read metadata files for visualization.")
    parser.add_argument("--input-json", required=True, help="Path to the *initial* list-based JSON file with anchor and read data.")
    parser.add_argument("--input-tsv", required=True, help="Path to the TSV file mapping anchor IDs to snarl IDs (must contain a header).")
    parser.add_argument("--output-json", required=True, help="Path for the output JSON file with transformed anchor data.")
    parser.add_argument("--snarl-output", required=True, help="Path for the output JSON file mapping snarls to anchors.")
    parser.add_argument("--connectivity-output", required=True, help="Path for the output JSON file showing anchor connectivity.")
    parser.add_argument("--read-metadata-output", required=True, help="Path for the output JSON file with metadata for each read.")
    args = parser.parse_args()

    try:
        with open(args.input_json, 'r') as f:
            raw_data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error reading input JSON file '{args.input_json}': {e}", file=sys.stderr)
        sys.exit(1)

    # Step 1: Transform the raw data and save to output JSON
    transformed_data = transform_raw_json(raw_data)
    with open(args.output_json, 'w') as f:
        json.dump(transformed_data, f, indent=4)
    print(f"Successfully transformed raw data and saved to '{args.output_json}'")

    # Step 2: Generate and save snarl-to-anchor mapping
    snarl_data = create_snarl_to_anchor_mapping(transformed_data, args.input_tsv)
    with open(args.snarl_output, 'w') as f:
        json.dump(snarl_data, f, indent=4)
    print(f"Successfully generated snarl mapping at '{args.snarl_output}'")

    # Step 3: Generate and save anchor connectivity data
    connectivity_data = create_anchor_connectivity(transformed_data)
    with open(args.connectivity_output, 'w') as f:
        json.dump(connectivity_data, f, indent=4)
    print(f"Successfully generated anchor connectivity data at '{args.connectivity_output}'")

    # Step 4: Generate and save read metadata
    read_metadata = create_read_metadata(transformed_data)
    with open(args.read_metadata_output, 'w') as f:
        json.dump(read_metadata, f, indent=4)
    print(f"Successfully generated read metadata at '{args.read_metadata_output}'")


if __name__ == "__main__":
    main()
