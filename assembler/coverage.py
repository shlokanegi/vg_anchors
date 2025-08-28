from bdsg.bdsg import PackedGraph

from assembler.gaf_reader import GafReader
import assembler.parser as lp

import time
import sys
from collections import defaultdict

from assembler.config import settings
from assembler.aligner import AlignAnchor


class Coverage(object):
    def __init__(self, packed_graph_path, gaf_path, min_cov=5) -> None:
        # useful initialization objects
        self.graph = PackedGraph()
        self.graph.deserialize(packed_graph_path)
        self.gaf_reader = GafReader(gaf_path)

        # important generated_data
        self.node_count =  defaultdict(int)
        self.min_coverage = min_cov
        self.aligner = AlignAnchor()
        self.aligner.build(dict_path, packed_graph_path)
    
    def get_total_basepairs(self, nodes):
        return sum(self.graph.get_length(self.graph.get_handle(node)) for node in nodes)
    
    def get_frequent_nodes(self, min_frequency):
        for node, count in self.node_count.items():
            if count >= min_frequency:
                yield node


    def get_alignment_coverage(self):
        for line in self.gaf_reader.get_lines():
            parsed_data = lp.processGafLine(line)
            if parsed_data:
                for node in parsed_data[settings.getint('NODE_POSITION')]:
                    self.node_count[node] += 1

            tot_bp = self.get_total_basepairs([x for x in self.get_frequent_nodes(self.min_coverage)])
        return tot_bp


if __name__ == "__main__":
    graph_path = sys.argv[1]
    alignment_path = sys.argv[2]
    coverage_obj = Coverage(graph_path, alignment_path)

    total_covered_bp = coverage_obj.get_alignment_coverage()

    print(f"Total base pairs covered by > {coverage_obj.min_coverage} : {total_covered_bp}")
