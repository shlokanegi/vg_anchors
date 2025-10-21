#!/usr/bin/env python3
"""
Script to parse snarl_tree_map.json and list all children (snarls and chains) of a given
parent ID along with their total leaf snarl counts.

Usage:
    python3 list_child_chains.py [input_file] [parent_id]
    
    input_file: Path to snarl_tree_map.json (default: snarl_tree_map.json)
    parent_id: The ID of the parent snarl/chain to query (default: root)
    
Examples:
    python3 list_child_chains.py                           # List root's children
    python3 list_child_chains.py snarl_tree_map.json root  # Same as above
    python3 list_child_chains.py snarl_tree_map.json c94412226-76657917  # List specific chain's children
"""

import json
import sys
from pathlib import Path

def main():
    # Default input file and parent ID
    input_file = "snarl_tree_map.json"
    parent_id = "root"
    
    # Parse command line arguments
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    if len(sys.argv) > 2:
        parent_id = sys.argv[2]
    
    if not Path(input_file).exists():
        print(f"Error: File '{input_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    print(f"Loading {input_file}...", file=sys.stderr)
    
    # Load the JSON file
    with open(input_file, 'r') as f:
        snarl_tree_map = json.load(f)
    
    # Get the parent entry
    if parent_id not in snarl_tree_map:
        print(f"Error: ID '{parent_id}' not found in the snarl tree map", file=sys.stderr)
        sys.exit(1)
    
    parent_children, parent_leaf_count = snarl_tree_map[parent_id]
    
    # Get all children (both snarls and chains)
    all_children = parent_children
    chains = [child for child in all_children if child.startswith('c')]
    snarls = [child for child in all_children if child.startswith('s')]
    
    print(f"\nQuerying parent: {parent_id}", file=sys.stderr)
    print(f"Found {len(all_children)} total children ({len(chains)} chains, {len(snarls)} snarls)", file=sys.stderr)
    print(f"Total leaf snarls at parent: {parent_leaf_count}\n", file=sys.stderr)
    
    # Create a list of (child_id, leaf_snarls, type) tuples
    child_data = []
    total_leaf_snarls = 0
    
    for child_id in all_children:
        if child_id not in snarl_tree_map:
            print(f"Warning: Child '{child_id}' not found in map", file=sys.stderr)
            continue
        
        children, leaf_snarls = snarl_tree_map[child_id]
        child_type = "chain" if child_id.startswith('c') else "snarl"
        child_data.append((child_id, leaf_snarls, child_type))
        total_leaf_snarls += leaf_snarls
    
    # Sort by leaf snarls (descending)
    child_data.sort(key=lambda x: x[1], reverse=True)
    
    # Print header
    print(f"{'ID':<35} {'Type':<8} {'Leaf Snarls':>15}")
    print("-" * 60)
    
    # Print each child
    for child_id, leaf_snarls, child_type in child_data:
        print(f"{child_id:<35} {child_type:<8} {leaf_snarls:>15,}")
    
    # Print summary
    print("-" * 60)
    print(f"{'TOTAL':<35} {'':<8} {total_leaf_snarls:>15,}")
    print(f"\nTotal children: {len(child_data)} ({len(chains)} chains, {len(snarls)} snarls)")

if __name__ == "__main__":
    main()

