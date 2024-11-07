from sys import argv

with open(argv[1], "r") as f:
    for line in f:
        clean_l = line.strip().split()
        name = clean_l[0]
        strand = clean_l[4]
        seq_start = int(clean_l[2])
        seq_end = int(clean_l[3])
        path_start = int(clean_l[7])
        path_end = int(clean_l[8])

        if strand == "-":
            print(f"read {name} has {strand} strand")
        if seq_end < seq_start:
            print(f"read {name} has seq_end < seq_start")
        if path_end < path_start:
            print(f"read {name} has path_end < path_start")