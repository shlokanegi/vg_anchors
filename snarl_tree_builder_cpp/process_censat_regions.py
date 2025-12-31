#!/usr/bin/env python3
"""
Process centromere/satellite BED file to identify coarse boundaries.
Step 1: Read CHM13 censat BED file and identify boundary nodes defining big enough centromere regions.

This script:
1. Reads the BED file and groups regions by chromosome
2. Calculates gaps between consecutive regions per chromosome
3. Merges close-by regions based on gap analysis
4. Outputs gap statistics and merged regions

Usage:
    python3 process_censat_regions.py [--input-bed PATH] [--output-gap-report PATH] [--output-merged-regions PATH]

Options:
    --input-bed PATH: Path to input BED file (default: chm13v2.0.cenSat.v2.0.bed)
    --output-gap-report PATH: Path to output gap report file (default: censat_gaps_report.txt)
    --output-merged-regions PATH: Path to output merged regions file (default: censat_merged_regions.bed)
"""

import sys
import argparse
from collections import defaultdict
import statistics

def parse_bed_file(bed_file):
    """
    Parse BED file and extract chromosome, start, end for each region.
    Returns a dictionary: {chromosome: [(start, end, annotation), ...]}
    """
    regions_by_chr = defaultdict(list)
    
    with open(bed_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip header line
            if line.startswith('track'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 4:
                continue
            
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            annotation = fields[3]
            
            regions_by_chr[chrom].append((start, end, annotation))
    
    # Sort regions by start position for each chromosome
    for chrom in regions_by_chr:
        regions_by_chr[chrom].sort(key=lambda x: x[0])
    
    return regions_by_chr

def calculate_gaps(regions_by_chr):
    """
    Calculate gaps between consecutive regions for each chromosome.
    Returns: {chromosome: [(gap_size, region1_end, region2_start), ...]}
    Also returns overall gap statistics.
    """
    gaps_by_chr = defaultdict(list)
    all_gaps = []
    
    for chrom, regions in sorted(regions_by_chr.items()):
        for i in range(len(regions) - 1):
            current_end = regions[i][1]
            next_start = regions[i+1][0]
            
            # Gap is the distance between end of current and start of next
            gap = next_start - current_end
            
            # Only consider positive gaps (non-overlapping regions)
            if gap >= 0:
                gaps_by_chr[chrom].append((gap, current_end, next_start, i, i+1))
                all_gaps.append(gap)
    
    return gaps_by_chr, all_gaps

def compute_gap_statistics(all_gaps):
    """
    Compute statistics on gap sizes to help determine merge threshold.
    """
    if not all_gaps:
        return {}
    
    sorted_gaps = sorted(all_gaps)
    stats = {
        'count': len(all_gaps),
        'min': min(all_gaps),
        'max': max(all_gaps),
        'mean': statistics.mean(all_gaps),
        'median': statistics.median(all_gaps),
        'stdev': statistics.stdev(all_gaps) if len(all_gaps) > 1 else 0,
        'percentile_10': sorted_gaps[int(len(sorted_gaps) * 0.10)],
        'percentile_25': sorted_gaps[int(len(sorted_gaps) * 0.25)],
        'percentile_50': sorted_gaps[int(len(sorted_gaps) * 0.50)],
        'percentile_75': sorted_gaps[int(len(sorted_gaps) * 0.75)],
        'percentile_90': sorted_gaps[int(len(sorted_gaps) * 0.90)],
        'percentile_95': sorted_gaps[int(len(sorted_gaps) * 0.95)],
        'percentile_99': sorted_gaps[int(len(sorted_gaps) * 0.99)],
    }
    
    return stats

def determine_merge_threshold(gap_stats):
    """
    Determine merge threshold based on gap statistics.
    
    Heuristic: We want to merge regions that are very close together (likely part of
    the same functional centromere) but keep separate regions that are far apart.
    
    Strategy:
    - Look at the distribution of gaps
    - Use a threshold around the 25th percentile or a fixed value (whichever is smaller)
    - This ensures we merge the closest regions while preserving distinct centromere blocks
    """
    if not gap_stats:
        return 100000  # Default 100kb
    
    # Use 25th percentile as a data-driven threshold
    # But cap it at 100kb to avoid merging truly separate regions
    threshold = min(gap_stats['percentile_25'], 100000)
    
    # If the 25th percentile is very small (< 10kb), use 50kb as a reasonable default
    if threshold < 10000:
        threshold = 50000
    
    return int(threshold)

def merge_regions(regions_by_chr, gaps_by_chr, merge_threshold):
    """
    Merge regions that are within merge_threshold distance of each other.
    Returns merged regions: {chromosome: [(start, end, region_count, annotations), ...]}
    """
    merged_by_chr = {}
    
    for chrom, regions in sorted(regions_by_chr.items()):
        if not regions:
            continue
        
        merged = []
        current_start = regions[0][0]
        current_end = regions[0][1]
        current_annotations = [regions[0][2]]
        region_count = 1
        
        for i in range(1, len(regions)):
            gap = regions[i][0] - current_end
            
            if gap <= merge_threshold:
                # Merge: extend the current region
                current_end = regions[i][1]
                current_annotations.append(regions[i][2])
                region_count += 1
            else:
                # Save current merged region and start a new one
                merged.append((
                    current_start,
                    current_end,
                    region_count,
                    ';'.join(set(current_annotations))  # Unique annotations
                ))
                
                current_start = regions[i][0]
                current_end = regions[i][1]
                current_annotations = [regions[i][2]]
                region_count = 1
        
        # Don't forget the last region
        merged.append((
            current_start,
            current_end,
            region_count,
            ';'.join(set(current_annotations))
        ))
        
        merged_by_chr[chrom] = merged
    
    return merged_by_chr

def write_gap_report(gaps_by_chr, all_gaps, gap_stats, output_file):
    """
    Write detailed gap analysis to file.
    """
    with open(output_file, 'w') as f:
        f.write("# Gap Analysis Report\n")
        f.write("# Analysis of distances between consecutive centromere/satellite regions\n\n")
        
        f.write("## Overall Gap Statistics\n")
        f.write(f"Total gaps: {gap_stats['count']}\n")
        f.write(f"Min gap: {gap_stats['min']:,} bp\n")
        f.write(f"Max gap: {gap_stats['max']:,} bp\n")
        f.write(f"Mean gap: {gap_stats['mean']:,.2f} bp\n")
        f.write(f"Median gap: {gap_stats['median']:,} bp\n")
        f.write(f"Std dev: {gap_stats['stdev']:,.2f} bp\n")
        f.write(f"10th percentile: {gap_stats['percentile_10']:,} bp\n")
        f.write(f"25th percentile: {gap_stats['percentile_25']:,} bp\n")
        f.write(f"50th percentile: {gap_stats['percentile_50']:,} bp\n")
        f.write(f"75th percentile: {gap_stats['percentile_75']:,} bp\n")
        f.write(f"90th percentile: {gap_stats['percentile_90']:,} bp\n")
        f.write(f"95th percentile: {gap_stats['percentile_95']:,} bp\n")
        f.write(f"99th percentile: {gap_stats['percentile_99']:,} bp\n\n")
        
        f.write("## Gaps by Chromosome\n")
        f.write("# Format: chromosome\tgap_size\tregion1_end\tregion2_start\tregion1_idx\tregion2_idx\n")
        
        for chrom, gaps in sorted(gaps_by_chr.items()):
            for gap, end1, start2, idx1, idx2 in gaps:
                f.write(f"{chrom}\t{gap}\t{end1}\t{start2}\t{idx1}\t{idx2}\n")

def write_merged_regions(merged_by_chr, output_file, merge_threshold):
    """
    Write merged centromere regions to BED-like file.
    """
    with open(output_file, 'w') as f:
        f.write(f"# Merged Centromere Regions (merge threshold: {merge_threshold:,} bp)\n")
        f.write("# Format: chromosome\tstart\tend\tregion_size\tnum_original_regions\tannotations\n")
        
        for chrom, regions in sorted(merged_by_chr.items()):
            for start, end, count, annotations in regions:
                size = end - start
                f.write(f"{chrom}\t{start}\t{end}\t{size}\t{count}\t{annotations}\n")

def main():
    # Currently hardcoded input and output file paths. TODO: make these command line arguments.
    input_bed = "chm13v2.0.cenSat.v2.0.bed"
    gap_report_file = "censat_gaps_report.txt"
    merged_regions_file = "censat_merged_regions.bed"
    
    print("Step 1: Parsing BED file...")
    regions_by_chr = parse_bed_file(input_bed)
    total_regions = sum(len(regions) for regions in regions_by_chr.values())
    print(f"  Found {total_regions} regions across {len(regions_by_chr)} chromosomes")
    
    print("\nStep 2: Calculating gaps between consecutive regions...")
    gaps_by_chr, all_gaps = calculate_gaps(regions_by_chr)
    print(f"  Calculated {len(all_gaps)} gaps")
    
    print("\nStep 3: Computing gap statistics...")
    gap_stats = compute_gap_statistics(all_gaps)
    print(f"  Median gap: {gap_stats['median']:,} bp")
    print(f"  Mean gap: {gap_stats['mean']:,.0f} bp")
    print(f"  25th percentile: {gap_stats['percentile_25']:,} bp")
    print(f"  75th percentile: {gap_stats['percentile_75']:,} bp")
    
    print("\nStep 4: Determining merge threshold...")
    merge_threshold = determine_merge_threshold(gap_stats)
    print(f"  Merge threshold: {merge_threshold:,} bp")
    print(f"  Rationale: Using min(25th percentile, 100kb) with 50kb minimum")
    print(f"            This merges the closest ~25% of region pairs")
    
    print("\nStep 5: Merging close-by regions...")
    merged_by_chr = merge_regions(regions_by_chr, gaps_by_chr, merge_threshold)
    total_merged = sum(len(regions) for regions in merged_by_chr.values())
    print(f"  Result: {total_merged} merged regions (from {total_regions} original)")
    print(f"  Reduction: {total_regions - total_merged} regions merged")
    
    print(f"\nStep 6: Writing gap report to {gap_report_file}...")
    write_gap_report(gaps_by_chr, all_gaps, gap_stats, gap_report_file)
    
    print(f"\nStep 7: Writing merged regions to {merged_regions_file}...")
    write_merged_regions(merged_by_chr, merged_regions_file, merge_threshold)
    
    print("\n" + "="*70)
    print("Summary:")
    print(f"  Input regions: {total_regions}")
    print(f"  Merged regions: {total_merged}")
    print(f"  Merge threshold: {merge_threshold:,} bp")
    print(f"  Gap report: {gap_report_file}")
    print(f"  Merged regions: {merged_regions_file}")
    print("="*70)

if __name__ == "__main__":
    main()




