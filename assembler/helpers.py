from sys import argv, stderr, exit
import json
from collections import defaultdict
from assembler.anchor import Anchor
from assembler.constants import RANGES, NUM_BINS, MIN_ANCHOR_LENGTH
import matplotlib.pyplot as plt
import pickle
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
        print(f"Error opening file {filename}: {e}")
        raise

def fastq_lines(in_fastqs):
    for fname in in_fastqs:
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


def plot_count_histogram(anchors_dict_fname: str, out_png: str) -> None:

    with open(anchors_dict_fname, 'rb') as in_f:
        sentinel_to_anchor = pickle.load(in_f)

    reads_count = defaultdict(int)
    for sentinel in sentinel_to_anchor:
        for anchor in sentinel_to_anchor[sentinel]:
            if anchor.num_sequences > 0:
                reads_count[anchor.num_sequences] += 1

    plt.bar(reads_count.keys(), reads_count.values())
    plt.xlabel("# Reads in anchors")
    plt.ylabel("Count")
    plt.title("# Reads in anchors distribution")
    plt.tight_layout()
    plt.savefig(out_png)


def plot_anchor_count_genome_distribution(anchors_dict_fname: str, out_png: str, title: str) -> None:

    count_dict = defaultdict(list)
    
    with open(anchors_dict_fname, 'rb') as in_f:
        sentinel_to_anchor = pickle.load(in_f)

    for sentinel in sentinel_to_anchor:
        for anchor in sentinel_to_anchor[sentinel]:
            position = anchor.genomic_position
            if position <= 0:
                continue
            count_dict[anchor.num_sequences].append(position)

    sorted_counts = sorted(count_dict.keys())
    print(f"{sorted_counts!r}")
    positions = [count_dict[count] for count in sorted_counts]

    print(positions)

    # Create the figure and axes
    fig, ax = plt.subplots(figsize=(24, 12))

    # Plot the stacked histogram
    ax.hist(positions, bins=NUM_BINS, stacked=True, label=sorted_counts)

    # Set title and labels
    ax.set_title(
        f'Anchor (size >={MIN_ANCHOR_LENGTH}) Count Distribution Across on {title}'
    )
    ax.set_xlabel(f"{title}")
    ax.set_ylabel("Number of Anchors")
    ax.legend(title="Reads count", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    binned_positions = []
    label = ["0"]
    binned_positions.append(count_dict[0])
    start_range = RANGES[0]

    for end_range in RANGES[1:]:
        new_range = []
        for count in range(start_range, end_range):
            if count_dict.get(count):
                new_range.extend(count_dict.get(count))
        binned_positions.append(new_range)
        label.append(f"[{start_range},{end_range})")
        start_range = end_range

    fig, ax = plt.subplots(figsize=(24, 12))

    # Plot the stacked histogram
    ax.hist(binned_positions, bins=NUM_BINS, stacked=True, label=label)

    # Set title and labels
    ax.set_title(
        f'Anchor (size >={MIN_ANCHOR_LENGTH}) Count Distribution Across {title}'
    )
    ax.set_xlabel(f"{title}")
    ax.set_ylabel("Number of Anchors")
    ax.legend(title="Reads count", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    plt.savefig(out_png[:-4] + ".binned.png", dpi=300, bbox_inches="tight")
    plt.close(fig)



def plot_heteroxigosity_on_genome(anchors_dict_fname: str, out_png: str, title: str) -> None:
    # import pkl sentinel_to_anchor_dictionary
    with open(anchors_dict_fname, 'rb') as in_f:
        sentinel_to_anchor = pickle.load(in_f)

    # import jsonl anchors
    # with open(anchors_json, "r") as f:
    #     anchors_file = json.load(f)

    ### 1 ###
    # scan the dictionary, if you find an anchor with > 1 read
    # populate a dictionary with key the snarl_id and item a list of tuples ("anchor_name", num_reads, position)
    # delete keys for snarl_id of just 1 tuple (genome not found heterozygous there)
    # now for every key check that the position of the anchors is the same, else take 1st
    heteroxygous_anchors = defaultdict(list)
    
    for sentinel in sentinel_to_anchor:
        for anchor in sentinel_to_anchor[sentinel]:
            if anchor.num_sequences >= 1:
                heteroxygous_anchors[anchor.snarl_id].append((repr(anchor),anchor.num_sequences, anchor.genomic_position))

    #removing omozygous loci
    # snarl_ids_to_remove = []
    positions_het = []
    positions_homo = []
    counts_het = []
    counts_homo = []
    for snarl_id, snarl_anchors in heteroxygous_anchors.items():
        print(f"visiting snarl {snarl_id}",end="\t")
        counts = []
        if len(snarl_anchors) > 1:
            snarl_positions = [x[2] for x in snarl_anchors]
            position = sum(snarl_positions) // len(snarl_positions)
            for element in snarl_anchors:
                positions_het.append(position)
                counts_het.append(element[1])
                counts.append(element[1])
            print(f"Position: {position}\tNum elements: { len(snarl_positions)}\tcounts: {counts!r}")
        else:
            positions_homo.append(snarl_anchors[0][2])
            counts_homo.append(snarl_anchors[0][1])
            print(f"Position: {snarl_anchors[0][2]}\tNum elements: 1\tcount: {snarl_anchors[0][1]}")


    # Create the figure and axes
    fig, ax = plt.subplots(figsize=(24, 12))

    # print(positions_homo)
    # print(counts_homo)
    # Plot the stacked histogram
    ax.scatter(positions_het, counts_het, label="Heterozygous")
    ax.scatter(positions_homo, counts_homo, label="Homozygous")

    # Set title and labels
    ax.set_title(
        f'Homo/hetero - zygous snarls for {title}'
    )
    ax.set_xlabel(f"Position in CHM13")
    ax.set_ylabel("Number of reads in each anchor")
    ax.legend(title="Zygosity") #, bbox_to_anchor=(1.05, 1), loc="upper left"
    plt.tight_layout()

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    ### 2 ###
    #plot a binned version of this graph just below

    ### 3 ###
    # for every key in the dict, create an array of the position and sort it.
    # violin plot of the distances between these values

    ### 4 ###
    # create dictionary of anchors to heterozygous or not (true or false)
    # create defaultdictionary with all the reads alinged to the graph
    # now for every time they are in an heterozygous , add +1 to the count of the defaultdict
    # plot the histogram on counts

    ### 5 ###
    # store all these info in a csv for a mutliplot on different values
    pass


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
