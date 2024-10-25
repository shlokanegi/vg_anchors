from bdsg.bdsg import HashGraph
from bdsg.bdsg import SnarlDistanceIndex
from bdsg.bdsg import PackedGraph

class SnarlAnchor:
    def __init__(self, max_nodes_in_snarl: int) -> None:
        self.graph = PackedGraph()
        self.index = SnarlDistanceIndex()
        self.leaf_snarls: list = []
        self.contains_child_snarls: bool = False
        self.num_nodes: int = 0
        self.max_num_nodes: int = max_nodes_in_snarl
        self.temp_nodes = []
        self.sentinel_to_anchor: dict = {}
    
    def build_graph(self, packed_graph_path: str,index_path: str )-> None:
        self.graph.deserialize(packed_graph_path)
        self.index.deserialize(index_path)

    def simple_snarl_iteratee(self,snarl_handle):
        print("Snarl:", self.index.net_handle_as_string(snarl_handle))
        return True

    def simple_chain_iteratee(self,chain_handle):
        print("Chain:", self.index.net_handle_as_string(chain_handle))
        return True

    def simple_node_iteratee(self,node_handle):
        print("Node:", self.index.net_handle_as_string(node_handle))
        return True
    
    def check_for_snarl(self,child_net_handle):

        if self.index.is_snarl(child_net_handle):
            self.contains_child_snarls = True
        elif self.index.is_node(child_net_handle):
            self.num_nodes += 1
            self.temp_nodes.append(child_net_handle)
        
        print(f"Children: {self.index.net_handle_as_string(child_net_handle)} has #nodes: {self.num_nodes} and has_snarls is: {self.contains_child_snarls!r}")
        return True
    
    def snarl_iteratee(self, net_handle):

        self.contains_child_snarls = False
        self.num_nodes = 0

        print(f"Visiting {self.index.net_handle_as_string(net_handle)}")

        snarl_children = []
        self.index.for_each_child(net_handle, lambda y: snarl_children.append(y) or True)

        self.temp_nodes = []
        for s_c in snarl_children:
            self.index.for_each_child(s_c, self.check_for_snarl)

        if not self.contains_child_snarls and self.num_nodes < self.max_num_nodes:
            self.leaf_snarls.append((net_handle, self.temp_nodes))

        return True
    
    def process_snarls(self):
        self.index.traverse_decomposition( 
            self.snarl_iteratee,    # snarl_iteratee
            lambda x: True, #  chain_iteratee
            lambda y: True  # node_iteratee
        )
        
        return self.leaf_snarls
    
    def print_tree_structure(self):
        self.index.traverse_decomposition( 
            self.simple_snarl_iteratee,    # snarl_iteratee
            self.simple_chain_iteratee, #  chain_iteratee
            self.simple_node_iteratee  # node_iteratee
        )