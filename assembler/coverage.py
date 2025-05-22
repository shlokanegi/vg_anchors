from bdsg.bdsg import PackedGraph

from assembler.gaf_reader import GafReader
import assembler.parser as lp

import time
import sys
from collections import defaultdict

from assembler.constants import (
    READ_P,
    R_LEN_P,
    STRAND_P,
    START_P,
    END_P,
    NODE_P,
    ORIENT_P,
    CS_P,
    READS_DEPTH,
)
            # read_name,
            # read_len,
            # relative_strand,
            # path_start,
            # path_end,
            # nodes_list,
            # orientation_list,
            # cs_line,

class GraphCoverage:
    def __init__(self, packed_graph_path, gaf_path, min_cov=5) -> None:
        # useful initialization objects
        self.graph = PackedGraph()
        self.graph.deserialize(packed_graph_path)
        self.gaf_reader = GafReader(gaf_path)

        # important generated_data
        self.node_count =  defaultdict(int)
        self.min_coverage = min_cov
    
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
                for node in parsed_data[NODE_P]:
                    self.node_count[node] += 1

            tot_bp = self.get_total_basepairs([x for x in self.get_frequent_nodes(self.min_coverage)])
        return tot_bp


if __name__ == "__main__":
    graph_path = sys.argv[1]
    alignment_path = sys.argv[2]
    coverage_obj = GraphCoverage(graph_path, alignment_path)

    total_covered_bp = coverage_obj.get_alignment_coverage()

    print(f"Total base pairs covered by > {coverage_obj.min_coverage} : {total_covered_bp}")
