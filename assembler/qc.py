from sys import argv, stderr
import json
import time
from assembler.rev_c import rev_c
import gzip
from contextlib import contextmanager
import itertools

READ_NAME_POS = 0
ORIENTATION_POS = 1
START_POS = 2
END_POS = 3


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
        print(fname)
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
                'header': header.strip(),
                'sequence': sequence.strip(),
                'plus_line': plus_line.strip(),
                'quality': quality.strip()
            }
        
        except StopIteration:
            break



def verify_anchors_validity(anchors_json: str, in_fastqs: list, out_fastq: str):

    with open(anchors_json, "r") as f:
        anchors_file = json.load(f)

    reads_ranges_dict = {}
    final_anchors_seq = []
    final_anchors_read_id = []

    ### FILL READS DICTIONARY WITH THE FOUND ANCHORS ###

    list_id = 0
    for anchor_list in anchors_file:
        for anchor in anchor_list[1:]:
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
        for _, start, end, anchor_id in values:
            ranges.append((start, end, anchor_id))
        soreted_ranges = sorted(ranges, key=lambda y: y[0])

        curr_start = -1
        curr_end = -1
        curr_id = 0
        for start, end, anchor_id in soreted_ranges:
            if start >= curr_end:
                curr_start = start
                curr_end = end
                curr_id = anchor_id
            else:
                print(
                    f"there is an overlap in read {read} \n between old anchor{curr_id} ({curr_start}-{curr_end}) and new anchor {anchor_id} ({start}-{end})",
                    file=stderr,
                )

    del anchor_list, soreted_ranges

    print(f"Found {len(final_anchors_seq)} sentinels used.", file=stderr)
    #print(f"Processing {in_fastq}", file=stderr)

    print_read = False
    count_fastq = 0
    read_count = 0

    ### VERIFY THAT ANCHORS IN THE SAME GROUP CORRESPOND TO THE SAME SUBSEQUENCES IN THE READS AND OUTPUT READS IN FASTQ ###
    for fname in in_fastqs:
        print(f"read: {fname}", file = stderr)

    with open(out_fastq, "w") as out_f:
        for entry in fastq_entries(fastq_lines(in_fastqs)):
            t0 = time.time()
            read_count += 1
            header = entry['header'].split('\t')[0] # MODIFIED FOR REVIO READS
            if reads_ranges_dict.get(header) != None:
                print_read = True
                #print(l, file=out_f)
                print(f"processing read {header}", end=" ", file=stderr)
                seq = entry['sequence']
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
            
        if print_read:
            print(f"{entry['header']}\n{entry['sequence']}\n{entry['plus_line']}\n{entry['quality']}", file=out_f)
            

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