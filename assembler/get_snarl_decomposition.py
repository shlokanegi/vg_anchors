from bdsg.bdsg import PackedGraph
from bdsg.bdsg import SnarlDistanceIndex
import time
import json


class SnarlTree:

    def __init__(self):
        self.graph = PackedGraph()
        self.index = SnarlDistanceIndex()
        self.snarl_tree_dict = {}   # the actual snarl tree dictionary(this is what we want to construct)
        self.net_to_position_in_snarl_tree = {}    # this is the path to the net in the snarl tree
        # self.net_to_name = {}    # this is the name of the net in the snarl tree
        self.current_cnt_of_snarl = 0
        self.current_cnt_of_chain = 0

    def build(self, graph_path, index_path, log_file):
        t0 = time.time()
        self.graph.deserialize(graph_path)
        print(f"Graph deserialized in {time.time()-t0:.2f}", file=log_file, flush=True)
        self.index.deserialize(index_path)
        print(f"Index deserialized in {time.time()-t0:.2f}", file=log_file, flush=True)

    def get_snarl_tree_dict(self, log_file):

        self.index.traverse_decomposition(
            snarl_iteratee=lambda s: self.func_snarl_iteratee(s, log_file),
            chain_iteratee=lambda c: self.func_chain_iteratee(c, log_file),
            node_iteratee=lambda n: self.func_node_iteratee(n, log_file)
        )
        return self.snarl_tree_dict
    
    def get_id(self, net):
        if self.index.is_chain(net) or self.index.is_snarl(net):
            left_bound_net = self.index.get_bound(net, False, False)
            right_bound_net = self.index.get_bound(net, True, False)
            left_bound_node_id = self.graph.get_id(self.index.get_handle(left_bound_net, self.graph))
            right_bound_node_id = self.graph.get_id(self.index.get_handle(right_bound_net, self.graph))
            if left_bound_node_id < right_bound_node_id:
                left_bound_node_id, right_bound_node_id = right_bound_node_id, left_bound_node_id
            bound_str = str(left_bound_node_id) + "-" + str(right_bound_node_id)
            if self.index.is_chain(net):
                return "c" + bound_str
            else:
                return "s" + bound_str
        elif self.index.is_root(net):
            return "root"
        else:
            raise ValueError(f"Net {net} is not a chain or snarl or root")
        
    def add_current_net_to_snarl_tree(self, net, current_net_id, log_file):
        parent_net = self.index.get_parent(net)
        parent_id = self.get_id(parent_net)
        if (parent_id not in self.net_to_position_in_snarl_tree):
            self.net_to_position_in_snarl_tree[parent_id] = "root"
            self.snarl_tree_dict["root"] = {}
        
        # print(f" ..Parent id: {parent_id}; Parent hash: {self.net_to_position_in_snarl_tree[parent_id]}", file=log_file, flush=True)
        
        root_to_net_path_in_snarl_tree = self.net_to_position_in_snarl_tree[parent_id].split("#") if self.net_to_position_in_snarl_tree[parent_id] is not None else []
        current_level_dict = self.snarl_tree_dict
        # looping to the parent dict (wiz. innermost dict) in the path
        for i in range(len(root_to_net_path_in_snarl_tree)):
            current_level_dict = current_level_dict[root_to_net_path_in_snarl_tree[i]]

        current_level_dict[current_net_id] = {}
        position_hash = ""
        if root_to_net_path_in_snarl_tree != []:
            position_hash = "#".join(root_to_net_path_in_snarl_tree) + "#" + str(current_net_id)
        else:
            position_hash = str(current_net_id)
        # print(f" ..Position hash: {position_hash}", file=log_file, flush=True)
        
        self.net_to_position_in_snarl_tree[current_net_id] = position_hash


    def func_node_iteratee(self, net, log_file):
        node_id = str(self.graph.get_id(self.index.get_handle(net,self.graph)))

        print(f"Node iteratee called with net {net} - node id {node_id}", file=log_file, flush=True)
        self.add_current_net_to_snarl_tree(net, node_id, log_file)
        return True

    def func_chain_iteratee(self, net, log_file):    
        chain_id = self.get_id(net)
        print(f"Chain iteratee called with net {net} - chain id {chain_id}", file=log_file, flush=True)
        self.add_current_net_to_snarl_tree(net, chain_id, log_file)
        return True
    
    def func_snarl_iteratee(self, net, log_file):
        snarl_id = self.get_id(net)
        print(f"Snarl iteratee called with net {net} - snarl id {snarl_id}", file=log_file, flush=True)
        self.add_current_net_to_snarl_tree(net, snarl_id, log_file)
        return True


def main():
    graph = "/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/results_hs/graph/PAW70337/PAW70337-16-sampled.pg.vg"
    # index = "/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/results_hs/hs-16/PAW70337/chr1_full/chunk/subgraph.dist"
    index = "/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/results_hs/hs-16/PAW70337/rccx/chunk/subgraph.pg.dist"
    
    # define outputs files (log file and snarl tree dict file)
    log_file_path = "snarl_tree_dict.log"
    snarl_tree_dict_file = "snarl_tree_dict.json"

    with open(log_file_path, "w") as log_file:
        snarl_tree = SnarlTree()
        snarl_tree.build(graph, index, log_file)
        snarl_tree_dict = snarl_tree.get_snarl_tree_dict(log_file)
    
    # save to json
    with open(snarl_tree_dict_file, "w") as f:
        json.dump(snarl_tree_dict, f, indent=4)

if __name__ == "__main__":
    main()
