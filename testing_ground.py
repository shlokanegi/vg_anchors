import sys
sys.path.append('./assembler')
#from anchor import SnarlAnchor
from path import parse_cs_line

# graph_path: str = 'small_test/chr20_small_idx.vg'
# index_path: str = 'small_test/chr20_small_idx.dist'

# anchoring = SnarlAnchor(10)
# anchoring.build_graph(graph_path, index_path)

# anchoring.print_tree_structure()

# leaf_snarl_net_handles: list = anchoring.get_leaf_snarls()

# one_snarl = leaf_snarl_net_handles[1]

# snarl_traversals = anchoring.get_paths_traversing_snarl(one_snarl)
l_count = 0

with open (sys.argv[1]) as f:
    for line in f:
        l = line.strip().split()
        if len(l) == 16:
            seq_len = int(l[1])
            equal = 0
            insertion = 0
            delition = 0
            subst = 0
            path_len = int(l[8]) - int(l[7])
            for el in parse_cs_line(l[14]):
                if el[0] == "+":
                    insertion += el[1]
                elif el[0] == ":":
                     equal += el[1]
                elif el[0] == "-":
                     delition += el[1]
                elif el[0] == "*":
                    subst += el[1]
                # print(f"F {el[0]} L {el[1]}", end=" ; ")
            tot_bases = equal + subst + delition

            print(f"At line {l_count} found alingment with estimated bases {tot_bases} but reported {path_len}.")
            print(f"Eq: {equal} ; Substitution: {subst} ; Insertions : {insertion} ; delitions : {delition}")
            print(f"Esitmated read length: {equal+subst+insertion} while reported {seq_len}")
            if l_count == 5:
                break
            l_count += 1
                            # elements.append(el[1])
            # print(f"\nAlignment Bases: {tot_bases-elements[0]}")
            # print(f"\nTotal Bases: {sum(tot_bases)}, alignment_bases:{sum(tot_bases)- tot_bases[0]}  ")