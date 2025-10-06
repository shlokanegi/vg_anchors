from sys import argv, stderr, exit
import json
from collections import defaultdict
from assembler.anchor import Anchor
from assembler.constants import MIN_ANCHOR_LENGTH, DEBUG
import gzip
from contextlib import contextmanager
from Bio import SeqIO

def reverse_complement(string) -> str:
    rev_str = string[::-1]
    r_c = ""
    for el in rev_str:
        if el == "A":
            r_c += "T"
        if el == "C":
            r_c += "G"
        if el == "G":
            r_c += "C"
        if el == "T":
            r_c += "A"
    return r_c

@contextmanager
def open_fastq(filename):
    try:
        if filename.endswith(".gz"):
            f = gzip.open(filename, 'rt')
        else:
            f = open(filename, 'r')
        try:
            yield f
        finally:
            f.close()
    except IOError as e:
        if DEBUG:
            print(f"Error opening file {filename}: {e}")
        raise

def fastq_lines(in_fastqs):
    for fname in in_fastqs:
        if DEBUG:
            print(fname,flush=True)
        with open_fastq(fname) as f:
            yield from f

def fastq_entries(fastq_lines_iter):
    """Generator that yields complete FASTQ entries"""
    while True:
        try:
            header = next(fastq_lines_iter)
            sequence = next(fastq_lines_iter)
            plus_line = next(fastq_lines_iter)
            quality = next(fastq_lines_iter)
            
            yield {
                'header': header.strip().split('\t')[0][1:],
                'sequence': sequence.strip(),
                'plus_line': plus_line.strip(),
                'quality': quality.strip()
            }
        
        except StopIteration:
            break


# Function to get complement
def complement(seq):
    # Define complement dictionary
    complement_map = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement_map)

# Function to get reverse complement
def rev_c(seq):
    return complement(seq)[::-1]

def extract_sequence(fasta_file, read_id):
    """
    get sequence for read from fasta
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == read_id:
            return str(record.seq)
    return None  # Return None if read_id is not found

if __name__ == "__main__":
    # verify_anchors_validity(argv[1], argv[2], argv[3])
    #anchors_shasta = argv[1]
    anchors_pos_dict = argv[1]
    out_png = argv[2]
    title = argv[3]

    plot_count_histogram(anchors_pos_dict, out_png + "count.png")

    plot_anchor_count_genome_distribution(anchors_pos_dict, out_png + "position_count.png", title
    )
