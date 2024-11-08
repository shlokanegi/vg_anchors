from sys import argv, stderr
import json
from collections import defaultdict
import matplotlib.pyplot as plt
import time
from assembler.rev_c import rev_c

NUM_BINS = 40



def plot_count_histogram(anchors_json: str, out_png: str) -> None:

    with open(anchors_json, "r") as f:
        anchors_list = json.load(f)
    
    reads_count = defaultdict(int)
    for anchor in anchors_list:
        reads_count[len(anchor)] += 1

    plt.bar(reads_count.keys(), reads_count.values())
    plt.xlabel("Reads in anchors")
    plt.ylabel("Count")
    plt.title(" Reads in anchors Frequency.")
    plt.tight_layout()
    plt.savefig(out_png)
    

def plot_anchor_count_genome_distribution(anchors_json: str, out_png: str) -> None: 
    # min_pos: int = -1
    # max_pos: int = -1
    count_dict = dict()

    with open(anchors_json, "r") as f:
        anchors_dict = json.load(f)

    for _ , anchor_l in anchors_dict.items():
        for anchor in anchor_l:
            (_, position, count) = anchor
            if position == -1:
                continue
            if count not in count_dict:
                count_dict[count] = []
            count_dict[count].append(position)
    
    sorted_counts = sorted(count_dict.keys())
    print(f"{sorted_counts!r}")
    positions = [count_dict[count] for count in sorted_counts]

    # Create the figure and axes
    fig, ax = plt.subplots(figsize=(24, 12))

    # Plot the stacked histogram
    ax.hist( positions, bins=NUM_BINS, stacked=True, label=sorted_counts)

    # Set title and labels
    ax.set_title('Anchor (size >=100) Count Distribution Across "CHM13#chr20:149948-250000')
    ax.set_xlabel('CHM13 Position 149948-250000')
    ax.set_ylabel('Number of Anchors')
    ax.legend(title='Reads count', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)

    binned_positions = []
    
    binned_positions.append( count_dict[0] )

    binned_positions.append([])
    for count in range(1,5):
        if count_dict.get(count):
            binned_positions[1].extend(count_dict.get(count))

    binned_positions.append([])
    for count in range(5,10):
        if count_dict.get(count):
            binned_positions[2].extend(count_dict.get(count))

    binned_positions.append([])
    for count in range(10,16):
        if count_dict.get(count):
            binned_positions[3].extend(count_dict.get(count))
    
    label = ['0', '[1,5)','[5,10)','[10,16)']

    fig, ax = plt.subplots(figsize=(24, 12))

    # Plot the stacked histogram
    ax.hist( binned_positions, bins=NUM_BINS, stacked=True, label=label)

    # Set title and labels
    ax.set_title('Anchor (size >=100) Count Distribution Across "CHM13#chr20:149948-250000')
    ax.set_xlabel('CHM13 Position in the interval 149948-250000')
    ax.set_ylabel('Number of Anchors')
    ax.legend(title='Reads count', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    plt.savefig(out_png[:-4]+'.binned.png', dpi=300, bbox_inches='tight')
    plt.close(fig)








def verify_anchors_validity(anchors_json: str, fastq: str, out_f: str):

    with open(anchors_json, "r") as f:
        anchors_list = json.load(f)

    reads_ranges_dict = {}
    final_anchors_seq = []
    final_anchors_read_id = []
    
    list_id = 0
    for anchor_l in anchors_list:
        for el in anchor_l:
            if reads_ranges_dict.get(el[0]) != None:
                reads_ranges_dict[el[0]].append(
                    [el[1], int(el[2]), int(el[3]), list_id]
                )
            else:
                reads_ranges_dict[el[0]] = [[el[1], int(el[2]), int(el[3]), list_id]]
        final_anchors_seq.append([])
        final_anchors_read_id.append([])
        list_id += 1

    for read, values in reads_ranges_dict.items():
        ranges = []
        for _, start, end, _ in values:
            ranges.append((start, end))
        soreted_ranges = sorted(ranges, key=lambda y: y[0])

        curr_start = -1
        curr_end = -1
        for start, end in soreted_ranges:
            if start >= curr_end:
                curr_start = start
                curr_end = end
            else:
                print(
                    f"there is an overlap in read {read} \n between old anchor {curr_start}-{curr_end} and new anchor {start}-{end}",
                    file=stderr,
                )

    del anchors_list
    read_c = 0
    print(f"Found {len(final_anchors_seq)} sentinels used.", file=stderr)
    print(f"Parsing {fastq}", file=stderr)
    print_read = False
    count_fastq = 0
    with open(fastq, "r") as f, open(out_f + ".fastq", "w") as out_fastq:
        for line in f:
            l = line.strip()
            t0 = time.time()
            if l[0] == "@":
                read_c += 1
                if reads_ranges_dict.get(l[1:]) != None:
                    print_read = True
                    print(l, file=out_fastq)
                    print(f"processing read {read_c}", end=" ", file=stderr)
                    sequence = f.readline()
                    seq = sequence.strip()
                    rev_s = rev_c(seq)
                    for elements in reads_ranges_dict.get(l[1:]):
                        if elements[1] >= len(seq) or elements[2] >= len(seq):
                            print(
                                f"Read {l[1:]} has len {len(seq)} but start {elements[1]} and end {elements[2]}"
                            )
                        if elements[0] == 1:
                            final_anchors_seq[elements[3]].append(
                                rev_s[elements[1] : elements[2]]
                            )
                        else:
                            final_anchors_seq[elements[3]].append(
                                seq[elements[1] : elements[2]]
                            )
                        final_anchors_read_id[elements[3]].append(
                            (l[1:], elements[0], elements[1], elements[2])
                        )
                    print(
                        f"in {time.time()-t0:.2f}. With {len(reads_ranges_dict.get(l[1:]))} elements.",
                        file=stderr,
                    )
                    print(seq, file=out_fastq)
            if print_read:
                if count_fastq == 3:
                    count_fastq = 0
                    print_read = False
                elif count_fastq > 0:
                    print(l, file=out_fastq)
                if print_read:
                    count_fastq += 1

    with open(out_f + ".id", "w") as outf:
        for line in final_anchors_read_id:
            print(f"{line!r}", file=outf)

    for idx, anchor in enumerate(final_anchors_seq):
        if len(anchor) > 0:
            if anchor.count(anchor[0]) != len(anchor):
                print(f"Anchors at {idx} do not match.", file=stderr)
            print(f"{anchor!r}")


if __name__ == "__main__":
    #verify_anchors_validity(argv[1], argv[2], argv[3])
    anchors_shasta = argv[1]
    anchors_count_pos_dict = argv[2]
    out_png = argv[3]

    plot_count_histogram(anchors_shasta, out_png+"count.png")

    plot_anchor_count_genome_distribution(anchors_count_pos_dict,out_png+"position_count.png")