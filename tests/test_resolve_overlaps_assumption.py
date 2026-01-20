#!/usr/bin/env python3
"""
Scan a GFA file and output link records (L lines) where:
1. A parent node has multiple children with *different* overlaps
2. A child node has multiple parents with *different* overlaps

Usage:
    python test.py -g input.gfa [--check parents|children|both]
    
Options:
    --check: Which checks to run (default: both)
        - parents: Check for parents with multiple children and different overlaps
        - children: Check for children with multiple parents and different overlaps
        - both: Run both checks (default)

If -g / --gfa is not provided, it defaults to the user's large GFA:
  /private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/\
results_hs/hs-16/PAW70337/wgtest_d50_D75_noCENT5000kb_80k_v1000_S10_a40/combined/\
combined_merged_assembly_p0.25_k16_c3.gfa
"""

import argparse
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-g",
        "--gfa",
        help="Input GFA file",
        default=(
            "/private/groups/migalab/shnegi/vg_anchors_project/"
            "test_lr_giraffe_assembly/results_hs/hs-16/PAW70337/"
            "wgtest_d50_D75_noCENT5000kb_80k_v1000_S10_a40/combined/"
            "combined_merged_assembly_p0.25_k16_c3.gfa"
        ),
    )
    parser.add_argument(
        "--check",
        choices=["parents", "children", "both"],
        default="both",
        help="Which checks to run (default: both)",
    )
    return parser.parse_args()


def find_parents_with_mixed_overlaps(gfa_path: str):
    """
    Read a GFA and print L-lines for parents that have multiple children
    with different overlaps.
    """
    # parent_id -> list of (line_str, overlap_str)
    parent_links = defaultdict(list)

    with open(gfa_path) as f:
        for line in f:
            if not line or line[0] != "L":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            parent = parts[1]
            overlap = parts[5]
            parent_links[parent].append((line.rstrip("\n"), overlap))

    found_any = False

    for parent, links in parent_links.items():
        overlaps = {ov for _, ov in links}
        # Only interested in parents with >= 2 children and >1 distinct overlap
        if len(links) >= 2 and len(overlaps) > 1:
            if not found_any:
                print("# Parents with multiple children and different overlaps")
            found_any = True
            print(f"# Parent: {parent} (overlaps: {', '.join(sorted(overlaps))})")
            for line, ov in links:
                print(line)
            print()

    if not found_any:
        print("# No parents found with multiple children and different overlaps.")


def find_children_with_mixed_overlaps(gfa_path: str):
    """
    Read a GFA and print L-lines for children that have multiple parents
    with different overlaps.
    """
    # child_id -> list of (line_str, overlap_str)
    child_links = defaultdict(list)

    with open(gfa_path) as f:
        for line in f:
            if not line or line[0] != "L":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            child = parts[3]  # to_id is the child
            overlap = parts[5]
            child_links[child].append((line.rstrip("\n"), overlap))

    found_any = False

    for child, links in child_links.items():
        overlaps = {ov for _, ov in links}
        # Only interested in children with >= 2 parents and >1 distinct overlap
        if len(links) >= 2 and len(overlaps) > 1:
            if not found_any:
                print("# Children with multiple parents and different overlaps")
            found_any = True
            print(f"# Child: {child} (overlaps: {', '.join(sorted(overlaps))})")
            for line, ov in links:
                print(line)
            print()

    if not found_any:
        print("# No children found with multiple parents and different overlaps.")


def main():
    args = parse_args()
    
    if args.check in ["parents", "both"]:
        find_parents_with_mixed_overlaps(args.gfa)
        if args.check == "both":
            print()
    
    if args.check in ["children", "both"]:
        find_children_with_mixed_overlaps(args.gfa)


if __name__ == "__main__":
    main()
