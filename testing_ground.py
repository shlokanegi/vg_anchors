import sys
sys.path.append('./assembler')
from flow_handler import Orchestrator
from anchor_dictionary_builder import SnarlAnchor
from line_parser import parse_cs_line
import time

graph_path: str = './large_test/chr20.full.100k.vg'
index_path: str = './large_test/chr20.full.100k.dist'
# dict_path: str = 'small_test/graph.dict'
sample_alignment: str = './large_test/m64012-190920-173625-Q20_chr20.full.100k.parsed.gaf'
anchors_file: str = './large_test/anchors_text.txt'
out_jsonl: str = './large_test/anchors_100k.json'
out_csv_bandage: str = './large_test/bandage_colors.csv'
out_d : str = './large_test/anchor_sizes.csv'

dictionary_builder = SnarlAnchor()
# Assume you have a method to build the dictionary
t0 = time.time()
dictionary_builder.build_graph(graph_path, index_path)
print(f"Indexes read in {time.time()-t0:.2f}", flush=True, file=sys.stderr)
# dictionary_builder.print_tree_structure()
t0 = time.time()
dictionary_builder.fill_anchor_sentinel_table()
print(f"Anchors dictionary built in {time.time()-t0:.2f}", flush=True, file=sys.stderr)
print(len(dictionary_builder.leaf_snarls))
dictionary = dictionary_builder.get_dict()
dictionary_builder.print_anchors_from_dict(anchors_file)
dictionary_builder.print_sentinels_for_bandage(out_csv_bandage)
dictionary_builder.print_dict_sizes(out_d)

orchestrator = Orchestrator(dictionary, graph_path, sample_alignment)

orchestrator.process()

orchestrator.dump_anchors(out_jsonl)

# t0 = time.time()
# time.sleep(5)
# print(f'Slept for {time.time()-t0:.2f} seconds.')
# one_snarl = leaf_snarl_net_handles[1]

#snarl_traversals = anchoring.get_paths_traversing_snarl(one_snarl)


# l_count = 0

# with open (sys.argv[1]) as f:
#     for line in f:
#         l = line.strip().split()
#         if len(l) == 16:
#             seq_len = int(l[1])
#             equal = 0
#             insertion = 0
#             delition = 0
#             subst = 0
#             path_len = int(l[8]) - int(l[7])
#             for el in parse_cs_line(l[14]):
#                 if el[0] == "+":
#                     insertion += el[1]
#                 elif el[0] == ":":
#                      equal += el[1]
#                 elif el[0] == "-":
#                      delition += el[1]
#                 elif el[0] == "*":
#                     subst += el[1]
#                 # print(f"F {el[0]} L {el[1]}", end=" ; ")
#             tot_bases = equal + subst + delition

#             print(f"At line {l_count} found alingment with estimated bases {tot_bases}; reported {path_len}.")
#             print(f"Eq: {equal} ; Substitution: {subst} ; Insertions : {insertion} ; delitions : {delition}")
#             print(f"Esitmated read length: {equal+subst+insertion} while reported {seq_len}")
#             if l_count == 5:
#                 break
#             l_count += 1
                            # elements.append(el[1])
            # print(f"\nAlignment Bases: {tot_bases-elements[0]}")
            # print(f"\nTotal Bases: {sum(tot_bases)}, alignment_bases:{sum(tot_bases)- tot_bases[0]}  ")

