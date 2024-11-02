import sys
sys.path.append('./assembler')
from flow_handler import Orchestrator
from anchor_dictionary_builder import SnarlAnchor
from line_parser import parse_cs_line

graph_path: str = './small_test/graph.vg'
index_path: str = './small_test/graph.dist'
dict_path: str = 'small_test/graph.dict'
sample_alignment: str = 'small_test/path_test.gaf'

dictionary_builder = SnarlAnchor()
# Assume you have a method to build the dictionary
dictionary_builder.build_graph(graph_path, index_path)
# dictionary_builder.print_tree_structure()
dictionary_builder.fill_anchor_sentinel_table()
print(len(dictionary_builder.leaf_snarls))
dictionary = dictionary_builder.get_dict()
dictionary_builder.print_anchors_from_dict()

orchestrator = Orchestrator(dictionary, graph_path, sample_alignment)

orchestrator.process()




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

