import sys
sys.path.append('./assembler')

from assembler.handler import Orchestrator
from assembler.builder import AnchorDictionary
import time
import os.path


def generate_dictionary(graph_path, index_path, anchors_dictionary):

        t0 = time.time()
        dictionary_builder = AnchorDictionary()
        dictionary_builder.build(graph_path, index_path)
        dictionary_builder.fill_anchor_dictionary()
        print(f"Anchors dictionary from {len(dictionary_builder.leaf_snarls)} snarls, containing {len(dictionary_builder.sentinel_to_anchor)} sentinels built in {time.time()-t0:.2f}", flush=True, file=sys.stderr)
        dictionary_builder.dump_dictionary(anchors_dictionary)

        # PRINT DEBUG INFO
        dictionary_builder.print_anchors_from_dict(anchors_file)
        dictionary_builder.print_sentinels_for_bandage(out_csv_bandage)
        dictionary_builder.print_dict_sizes(out_d)
        dictionary_builder.generate_positioned_dictionary('',out_positioned_dict)

