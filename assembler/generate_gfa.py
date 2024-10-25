#!/usr/bin/env python3.11
from sys import argv

def write_gfa_graph(file_out, nodes_seq, edges, paths):
    with open(file_out,"w") as f:
        count = 1
        for seq in nodes_seq:
            print(f'S\t{count}\t{seq}', file=f)
            count += 1
        for edge in edges:
            print(f'L\t{edge[0]}\t+\t{edge[1]}\t+\t0M', file=f)
        count = 1
        for path in paths:
            tot_l=0
            end_str=""
            for node in path:
                tot_l += len(nodes_seq[node-1])
                end_str += f">{node}"
            print(f'W\tP{count}\t{count}\tchr1\t0\t{tot_l}\t{end_str}', file=f)
            count += 1

if __name__ == "__main__":
    file_w = argv[1]
    nodes_seq = ["AAAAAAAA","BBBB","CCCC","DDDDD","EE","FFFF","GGG","H","I","L"]
    edges = [(1,2),(1,3),(3,4),(3,5),(4,7),(4,8),(7,5),(8,5),(5,9),(5,10),(9,6),(10,6),(6,2)]
    paths = [(1,3,5,10,6,2),(1,2),(1,2),(1,3,4,7,5,9,6,2),(1,3,4,7,5,9,6,2),(1,3,4,8,5,10,6,2),(1,3,4,8,5,10,6,2),(1,3,4,7,5,10,6,2)]

    write_gfa_graph(file_w, nodes_seq, edges, paths)


    # print(f'S\t1\tAAAA', file=f)
    # print(f'S\t2\tBB', file=f)
    # print(f'S\t3\tC', file=f)
    # print(f'S\t4\tDDDDD', file=f)
    # print(f'L\t1\t+\t2\t+\t0M', file=f)
    # print(f'L\t1\t+\t3\t+\t0M', file=f)
    # print(f'L\t2\t+\t4\t+\t0M', file=f)
    # print(f'L\t3\t+\t4\t+\t0M', file=f)
    # print(f'W\tP1\t1\tchr1\t0\t10\t>1>2>4', file=f)
    # print(f'W\tP2\t2\tchr1\t0\t9\t>1>3>4', file=f)

    # nodes_seq = ["AAAAAAAA","BBBB","CCCC","DDDDD","EE","FFFF","GGG","H","I"]
    # edges = [(1,2),(1,3),(3,6),(3,7),(6,4),(7,4),(4,8),(4,9),(8,5),(9,5),(5,2)]
    # paths = [(1,3,6,4,8,5,2),(1,2),(1,2),(1,3,6,4,8,5,2),(1,3,7,4,9,5,2),(1,3,7,4,9,5,2),(1,3,7,4,8,5,2)]

    # nodes_seq = ["AAAAAAAA","BBB","CCCC","DDDDD","EE","FFFF","GGG","H","I"]
    # edges = [(1,2),(2,3),(2,4),(3,5),(4,5),(5,6),(5,7),(7,8),(6,8),(8,9),(1,9)]
    # paths = [(1,2,3,5,7,8,9),(1,9),(1,9),(1,2,3,5,7,8,9),(1,2,4,5,6,8,9),(1,2,4,5,6,8,9),(1,2,3,4,5,8,9)]


    # nodes_seq = ["AAAAA","BB","C","DDDDD"]
    # edges = [(1,2),(1,3),(2,4),(3,4)]
    # paths = [(1,2,4),(1,2,4),(1,3,4)]