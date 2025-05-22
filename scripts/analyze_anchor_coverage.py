#!/usr/bin/env python3

import json
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Tuple
import argparse

def load_coverage_data(coverage_file: str) -> Tuple[Dict[str, int], Dict[str, int]]:
    """Load coverage statistics from JSON file"""
    with open(coverage_file, 'r') as f:
        data = json.load(f)
        return data['initial_coverage'], data['final_coverage']

def plot_coverage_histogram(initial: Dict[str, int], final: Dict[str, int], output_file: str):
    """Plot histograms of initial and final coverage"""
    plt.figure(figsize=(12, 6))
    
    # Plot initial coverage
    plt.subplot(1, 2, 1)
    plt.hist(list(initial.values()), bins=50, alpha=0.5, label='Initial Coverage')
    plt.title('Initial Anchor Coverage Distribution')
    plt.xlabel('Number of Reads')
    plt.ylabel('Number of Anchors')
    plt.legend()
    
    # Plot final coverage
    plt.subplot(1, 2, 2)
    plt.hist(list(final.values()), bins=50, alpha=0.5, label='Final Coverage')
    plt.title('Final Anchor Coverage Distribution')
    plt.xlabel('Number of Reads')
    plt.ylabel('Number of Anchors')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def print_coverage_stats(initial: Dict[str, int], final: Dict[str, int]):
    """Print coverage statistics"""
    print("\nCoverage Statistics:")
    print("-" * 50)
    
    # Initial coverage stats
    initial_values = list(initial.values())
    print("\nInitial Coverage:")
    print(f"Total anchors: {len(initial)}")
    print(f"Mean coverage: {np.mean(initial_values):.2f}")
    print(f"Median coverage: {np.median(initial_values):.2f}")
    print(f"Min coverage: {min(initial_values)}")
    print(f"Max coverage: {max(initial_values)}")
    
    # Final coverage stats
    final_values = list(final.values())
    print("\nFinal Coverage:")
    print(f"Total anchors: {len(final)}")
    print(f"Mean coverage: {np.mean(final_values):.2f}")
    print(f"Median coverage: {np.median(final_values):.2f}")
    print(f"Min coverage: {min(final_values)}")
    print(f"Max coverage: {max(final_values)}")
    
    # Coverage change stats
    common_anchors = set(initial.keys()) & set(final.keys())
    coverage_changes = [final[anchor] - initial[anchor] for anchor in common_anchors]
    print("\nCoverage Changes (Final - Initial):")
    print(f"Mean change: {np.mean(coverage_changes):.2f}")
    print(f"Median change: {np.median(coverage_changes):.2f}")
    print(f"Min change: {min(coverage_changes)}")
    print(f"Max change: {max(coverage_changes)}")

def main():
    parser = argparse.ArgumentParser(description='Analyze anchor coverage statistics')
    parser.add_argument('coverage_file', help='Path to coverage JSON file')
    parser.add_argument('--output', '-o', help='Output plot file path', default='coverage_histogram.png')
    args = parser.parse_args()
    
    # Load coverage data
    initial_coverage, final_coverage = load_coverage_data(args.coverage_file)
    
    # Print statistics
    print_coverage_stats(initial_coverage, final_coverage)
    
    # Plot histograms
    plot_coverage_histogram(initial_coverage, final_coverage, args.output)
    print(f"\nCoverage histogram saved to {args.output}")

if __name__ == '__main__':
    main() 
