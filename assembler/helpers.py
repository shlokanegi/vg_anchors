from sys import argv, stderr, exit
import json
from collections import defaultdict
from assembler.anchor import Anchor
from assembler.config import settings
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


def plot_count_histogram(pkl_file, out_png):
    import matplotlib.pyplot as plt
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)

    # data = sorted(data, key=lambda x: x[0])
    counts = [d[1] for d in data]

    plt.hist(counts, bins=settings.getint('NUM_BINS'), range=(min(counts), max(counts)), edgecolor="black")
    plt.xlabel("Count")
    plt.ylabel("Frequency")
    plt.title("Count distribution")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def plot_anchor_count_genome_distribution(pkl_file, out_png, title):
    import matplotlib.pyplot as plt
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)

    positions = [d[0] for d in data]
    counts = [d[1] for d in data]

    plt.figure(figsize=(20, 10))
    plt.bar(positions, counts, width=1.0)
    plt.xlabel("Position")
    plt.ylabel("Count")
    plt.title(f"Anchor count distribution across {title}")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()



def plot_heteroxigosity_on_genome(anchor_info_file, out_png, title):
    import matplotlib.pyplot as plt
    with open(anchor_info_file, "rb") as f:
        anchor_info = pickle.load(f)

    # import jsonl anchors
    # with open(anchors_json, "r") as f:
    #     anchors_file = json.load(f)

    ### 1 ###
    # scan the dictionary, if you find an anchor with > 1 read
    # populate a dictionary with key the snarl_id and item a list of tuples ("anchor_name", num_reads, position)
    # delete keys for snarl_id of just 1 tuple (genome not found heterozygous there)
    # now for every key check that the position of the anchors is the same, else take 1st
    heteroxygous_anchors = defaultdict(list)
    
    for sentinel in anchor_info:
        for anchor in anchor_info[sentinel]:
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
