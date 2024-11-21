import sys
import gzip
from contextlib import contextmanager

@contextmanager
def open_fasta(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, 'rt') as f:
            yield f
    else:
        with open(filename, 'r') as f:
            yield f

def read_bedfile(filename, target_chromosome):
    intervals = []
    with open(filename, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom, start, end = fields[:3]
            if chrom == target_chromosome:
                intervals.append((int(start), int(end)))
    return intervals

def merge_intervals(intervals):
    if not intervals:
        return []
    
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    
    for current in intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    
    return merged

def find_gaps(intervals, chromosome_length):
    gaps = []
    last_end = 0
    for start, end in intervals:
        if start > last_end:
            gaps.append((last_end, start))
        last_end = max(last_end, end)
    if last_end < chromosome_length:
        gaps.append((last_end, chromosome_length))
    return gaps

def get_size(genome_fasta, chromosome_name):
    chr_size = 0
    get_size = False
    with open_fasta(genome_fasta) as f:
        for line in f:
            if get_size and line.strip().startswith(">"):
                return chr_size
            elif get_size:
                chr_size += len(line.strip())
                continue
            if line.strip().startswith(f">{chromosome_name}"):
                get_size = True



if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script.py <bedfile1> <bedfile2> <bedfile3> <chromosome_name> <chromosome_length>")
        sys.exit(1)

    genome_fasta = sys.argv[1]
    bedfiles = sys.argv[2:-1]
    target_chromosome = sys.argv[-1]
    chromosome_length = get_size(genome_fasta, target_chromosome)

    print(f"Found {chromosome_length} as max length of {target_chromosome}.", file=sys.stderr)

    all_intervals = []
    for bedfile in bedfiles:
        all_intervals.extend(read_bedfile(bedfile, target_chromosome))

    merged_intervals = merge_intervals(all_intervals)
    gaps = find_gaps(merged_intervals, chromosome_length)

    # Print the gaps
    print(f"Found {len(gaps)} ranges not spanning the regions in the bedfiles for chromosome {target_chromosome}.", file=sys.stderr)
    for start, end in gaps:
        print(f"{target_chromosome}\t{start}\t{end}")