from sys import argv, stderr
import json
from collections import defaultdict
import matplotlib.pyplot as plt
import time

def verify_anchors_validity(anchors_json: str, fastq: str):

    with open(anchors_json, 'r') as f:
        anchors_list = json.load(f)
    
    reads_ranges_dict = {}
    final_anchors_seq = []
    anchor_multiplicity = defaultdict(int)
    list_id = 0
    for anchor_l in anchors_list:
        anchor_multiplicity[len(anchor_l)] += 1
        for el in anchor_l:
            if reads_ranges_dict.get(el[0]) != None:
                reads_ranges_dict[el[0]].append([int(el[2]),int(el[3]), list_id])
            else:
                reads_ranges_dict[el[0]] = [[int(el[2]),int(el[3]), list_id]]
        final_anchors_seq.append([])
        list_id += 1

    for read, values in reads_ranges_dict.items():
        ranges = []
        for start, end, _ in values:
            ranges.append((start, end))
        soreted_ranges = sorted(ranges, key = lambda y : y[0])

        curr_start = -1
        curr_end = -1
        for start, end in soreted_ranges:
            if start >= curr_end:
                curr_start = start
                curr_end = end
            else:
                print(f"there is an overlap in read {read} \n between old anchor {curr_start}-{curr_end} and new anchor {start}-{end}", file=stderr)




    plt.bar(anchor_multiplicity.keys(), anchor_multiplicity.values())
    plt.xlabel('Reads in anchors')
    plt.ylabel('Count')
    plt.title(' Reads in anchors Frequency.')
    plt.show()
    plt.savefig('./dict_plot.png')
    
    del anchors_list
    read_c = 0
    print(f"Found {len(final_anchors_seq)} sentinels used.", file=stderr)
    print(f"Parsing {fastq}", file=stderr)
    with open(fastq, "r") as f:
        for line in f:
            l = line.strip()
            t0 = time.time()
            if l[0] == "@":
                read_c += 1
                if reads_ranges_dict.get(l[1:]) != None:
                    print(f"processing read {read_c}", end=" ", file=stderr)
                    sequence = f.readline()
                    seq = sequence.strip()
                    for elements in reads_ranges_dict.get(l[1:]):
                        final_anchors_seq[elements[2]].append(seq[elements[0]:elements[1]])
                    print(f"in {time.time()-t0:.2f}. With {len(reads_ranges_dict.get(l[1:]))} elements.", file=stderr)
        
    
    for anchor in final_anchors_seq:
        if len(anchor) > 0:
            if anchor.count(anchor[0]) != len(anchor):
                print(f"Anchors do not match.", file=stderr)
            print(f"{anchor!r}")


if __name__ == "__main__":
    verify_anchors_validity(argv[1], argv[2])