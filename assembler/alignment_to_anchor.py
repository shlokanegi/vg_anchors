from bdsg.bdsg import PackedGraph

import pickle
from math import ceil, floor

class AlignAnchor:

    def __init__(self, packed_graph_path: str, anchor_dictionary_path: str) -> None:
        #useful initialization objects
        self.graph = PackedGraph()
        self.graph.deserialize(packed_graph_path)
        self.sentinel_to_anchor: dict = {}
        with open(anchor_dictionary_path, 'rb') as f:
            self.sentinel_to_anchor = pickle.load(f)

        # temporary variables to store data between functions

    