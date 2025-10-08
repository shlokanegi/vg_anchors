import argparse
import csv
import sys

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

def main():
    """Main function to run the analysis script."""
    parser = argparse.ArgumentParser(
        description="Outputs journey data for reads from a Journeys.csv file that contain anchors with both '+' and '-' orientations. The output is in CSV format."
    )
    parser.add_argument(
        "journeys_file", 
        help="Path to the Journeys.csv file."
    )
    args = parser.parse_args()

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

        writer = csv.writer(sys.stdout)
        writer.writerows(sorted_rows)
    else:
        print("No reads with mixed anchor strands were found.", file=sys.stderr)

if __name__ == "__main__":
    main()
