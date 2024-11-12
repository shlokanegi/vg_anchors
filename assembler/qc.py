from sys import argv, stderr
import json
import time
from assembler.rev_c import rev_c
import gzip
from contextlib import contextmanager

READ_NAME_POS = 0
ORIENTATION_POS = 1
START_POS = 2
END_POS = 3


@contextmanager
def open_fastq(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename) as f:
            yield f
    else:
        with open(filename) as f:
            yield f


def verify_anchors_validity(anchors_json: str, in_fastq: str, out_fastq: str):

    with open(anchors_json, "r") as f:
        anchors_file = json.load(f)

    reads_ranges_dict = {}
    final_anchors_seq = []
    final_anchors_read_id = []

    ### FILL READS DICTIONARY WITH THE FOUND ANCHORS ###

    list_id = 0
    for anchor_list in anchors_file:
        for anchor in anchor_list:
            sequence_name = anchor[READ_NAME_POS]
            orientation = int(anchor[ORIENTATION_POS])
            range_start = int(anchor[START_POS])
            range_end = int(anchor[END_POS])
            if reads_ranges_dict.get(sequence_name) != None:
                
                reads_ranges_dict[sequence_name].append(
                    [orientation, range_start, range_end, list_id]
                )
            else:
                reads_ranges_dict[sequence_name] = [[orientation, range_start, range_end, list_id]]
        final_anchors_seq.append([])
        final_anchors_read_id.append([])
        list_id += 1


    ### VERIFY THAT ANCHORS DO NOT OVERLAP IN THE SEQUENCE ###

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

    del anchor_list, soreted_ranges

    print(f"Found {len(final_anchors_seq)} sentinels used.", file=stderr)
    print(f"Processing {in_fastq}", file=stderr)

    print_read = False
    count_fastq = 0
    read_count = 0

    ### VERIFY THAT ANCHORS IN THE SAME GROUP CORRESPOND TO THE SAME SUBSEQUENCES IN THE READS AND OUTPUT READS IN FASTQ ###

    with open_fastq(in_fastq) as f, open(out_fastq, "w") as out_f:
        for line in f:
            l = line.strip()
            t0 = time.time()
            if l[0] == "@":
                read_count += 1
                header = l[1:].split('\t')[0] # MODIFIED FOR REVIO READS
                if reads_ranges_dict.get(header) != None:
                    print_read = True
                    print(l, file=out_f)
                    print(f"processing read {read_count}", end=" ", file=stderr)
                    sequence = f.readline()
                    seq = sequence.strip()
                    rev_s = rev_c(seq)
                    for elements in reads_ranges_dict.get(header):
                        if elements[1] >= len(seq) or elements[2] >= len(seq):
                            print(
                                f"Read {header} has len {len(seq)} but start {elements[1]} and end {elements[2]}", file=stderr
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
                            (header, elements[0], elements[1], elements[2])
                        )
                    print(
                        f"in {time.time()-t0:.2f}. With {len(reads_ranges_dict.get(header))} elements.",
                        file=stderr,
                    )
                    print(seq, file=out_f)
            if print_read:
                if count_fastq == 3:
                    count_fastq = 0
                    print_read = False
                elif count_fastq > 0:
                    print(l, file=out_f)
                if print_read:
                    count_fastq += 1

    with open(out_fastq + ".id", "w") as outf:
        for line in final_anchors_read_id:
            print(f"{line!r}", file=outf)

    for idx, anchor in enumerate(final_anchors_seq):
        if len(anchor) > 0:
            if anchor.count(anchor[0]) != len(anchor):
                print(f"Anchors at {idx} do not match.", file=stderr)
                print(f"{anchor!r}",file=stderr)

if __name__ == "__main__":
    anchors_shasta = argv[1]
    in_fastq = argv[2]
    out_fastq = in_fastq.rstrip(".fastq") + "selected.fastq" if in_fastq.endswith(".fastq") else in_fastq.rstrip(".fastq.gz") + "selected.fastq"

    verify_anchors_validity(anchors_shasta, in_fastq, out_fastq)