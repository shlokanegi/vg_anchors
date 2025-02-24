from sys import argv, stderr, exit
import json
from collections import defaultdict
from assembler.anchor import Anchor
from assembler.constants import RANGES, NUM_BINS, MIN_ANCHOR_LENGTH
import matplotlib.pyplot as plt
import pickle



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
        f'Anchor (size >=50) Count Distribution Across on {title}'
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
        binned_positions.append([])
        for count in range(start_range, end_range):
            if count_dict.get(count):
                binned_positions[1].extend(count_dict.get(count))
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
    positions = []
    counts = []
    for snarl_id, snarl_anchors in heteroxygous_anchors.items():
        print(f"visiting snarl {snarl_id}")
        if len(snarl_anchors) > 1:
            snarl_positions = [x[2] for x in snarl_anchors]
            position = sum(snarl_positions) // len(snarl_positions)
            for element in snarl_anchors:
                positions.append(position)
                counts.append(element[1])

    print(positions)
    print(counts)

    # Create the figure and axes
    fig, ax = plt.subplots(figsize=(24, 12))

    # Plot the stacked histogram
    ax.plot(positions, counts)

    # Set title and labels
    ax.set_title(
        f'Heteorozigous snarls for {title}'
    )
    ax.set_xlabel(f"Position in CHM13")
    ax.set_ylabel("Number of reads in each anchor")
    # ax.legend(title="Reads count", bbox_to_anchor=(1.05, 1), loc="upper left")
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


if __name__ == "__main__":
    # verify_anchors_validity(argv[1], argv[2], argv[3])
    #anchors_shasta = argv[1]
    anchors_pos_dict = argv[1]
    out_png = argv[2]
    title = argv[3]

    plot_count_histogram(anchors_pos_dict, out_png + "count.png")

    plot_anchor_count_genome_distribution(anchors_pos_dict, out_png + "position_count.png", title
    )
