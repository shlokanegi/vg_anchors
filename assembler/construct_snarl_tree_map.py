import os
import multiprocessing
from bdsg.bdsg import PackedGraph
from bdsg.bdsg import SnarlDistanceIndex
import time
import json
from sys import stderr

PARALLELIZATION_DEPTH = 2

class SnarlDecomposition:

    def __init__(self):
        self.graph = PackedGraph()
        self.index = SnarlDistanceIndex()
        self.snarl_tree_map = {}   # snarl_tree_map[snarl_id] = {[children_ids], count_of_leaf_snarls_in_subtree}

    def build(self, graph_path, index_path, log_file):
        t0 = time.time()
        self.graph.deserialize(graph_path)
        print(f"Graph deserialized in {time.time()-t0:.2f}", file=stderr, flush=True)
        self.index.deserialize(index_path)
        print(f"Index deserialized in {time.time()-t0:.2f}", file=stderr, flush=True)


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
        elif self.index.is_node(net):
            return str(self.graph.get_id(self.index.get_handle(net, self.graph)))
        else:
            raise ValueError(f"Net {net} is not a chain, snarl, node or root")

    def check_snarl_in_children_iteratee(self, child_net_handle) -> bool:

        if self.index.is_snarl(child_net_handle):
            self.contains_child_snarls = True
            return False
        return True

    def is_leaf_snarl(self, net):
        self.contains_child_snarls = False

        snarl_children: list = []
        self.index.for_each_child(
            net, lambda y: snarl_children.append(y) or True
        )

        for snarl_child in snarl_children:
            self.index.for_each_child(snarl_child, self.check_snarl_in_children_iteratee)

        if not self.contains_child_snarls:
            return True

    
    def traverse_decomposition_helper(self, handle, log_file) -> tuple[int, str]:
        map_key = self.get_id(handle)
        
        # The stopping condition is when the current net is a leaf snarl
        if self.index.is_snarl(handle) and self.is_leaf_snarl(handle):
            # self.snarl_tree_map[map_key] = ([], 1)
            return 1, map_key
        elif self.index.is_node(handle):
            return 0, map_key

        cnt_leaf_snarls_in_subtree = 0
        children_handles = []
        children_map_keys = []
        self.index.for_each_child(handle, lambda y: (children_handles.append(y) if not self.index.is_node(y) else True) or True)
        for child_handle in children_handles:
            cnt_leaf_snarls_in_child_subtree, child_map_key = self.traverse_decomposition_helper(child_handle, log_file)
            children_map_keys.append(child_map_key)
            # Once the recursive call for a child returns, it will have returned the total number of leaf snarls in that child's entire subtree. This count is added to the current snarl's running total.
            cnt_leaf_snarls_in_subtree += cnt_leaf_snarls_in_child_subtree
        
        # After the for loop ends, the function has visited all the current snarl's children and has received the total leaf snarl counts for each of them.
        self.snarl_tree_map[map_key] = (children_map_keys, cnt_leaf_snarls_in_subtree)
        
        return cnt_leaf_snarls_in_subtree, map_key


    def traverse_decomposition(self, subtree_root_handle, log_file):
        self.traverse_decomposition_helper(subtree_root_handle, log_file)


def main():
    graph = "/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/results_hs/hs-16/PAW70337/PDPK/query/subgraph.pg.vg"
    index = "/private/groups/migalab/shnegi/vg_anchors_project/test_lr_giraffe_assembly/results_hs/hs-16/PAW70337/PDPK/query/subgraph.pg.dist"
    
    # define outputs files (log file and snarl tree dict file)
    log_file_path = "snarl_tree_map.log"
    output_json_path = "snarl_tree_map.json"

    with open(log_file_path, "w") as log_file:
        obj = SnarlDecomposition()
        # Load the graph and index
        obj.build(graph, index, log_file)
        # Get the root handle
        root_handle = obj.index.get_root()
        # top_level_traversal_before_parallelization(root_handle)

        obj.traverse_decomposition(root_handle, log_file)


    # save to json
    with open(output_json_path, "w") as f:
        json.dump(obj.snarl_tree_map, f, indent=4)

if __name__ == "__main__":
    main()




# def traverse_decomposition_helper(handle):
#     map_key = get_id(handle)
    
#     # The stopping condition is when the current net is a leaf snarl
#     if_leaf_snarl = index.is_leaf_snarl(handle)
#     if is_leaf_snarl:    # base case
#         snarl_tree_map[map_key] = {[], 1}
#         return 1

#     cnt_leaf_snarls_in_subtree = 0
#     children_ids = []
#     for child_handle in index.get_children(handle):
#         child_id = get_id(child_handle)
#         children_ids.append(child_id)
#         cnt_leaf_snarls_in_child_subtree = traverse_decomposition_helper(child_handle)
#         # Once the recursive call for a child returns, it will have returned the total number of leaf snarls in that child's entire subtree. This count is added to the current snarl's running total.
#         cnt_leaf_snarls_in_subtree += cnt_leaf_snarls_in_child_subtree
    
#     # After the for loop ends, the function has visited all the current snarl's children and has received the total leaf snarl counts for each of them.
#     snarl_tree_map[map_key] = {children_ids, cnt_leaf_snarls_in_subtree}
#     return cnt_leaf_snarls_in_subtree


# def traverse_decomposition(subtree_root_handle):
#     global snarl_tree_map = {}    # snarl_tree_map[snarl_id] = {[children_ids], count_of_leaf_snarls_in_subtree}
#     traverse_decomposition_helper(subtree_root_handle)
#     dump_snarl_tree_map()


# def top_level_traversal_helper(handle, current_depth):
#     if current_depth == PARALLELIZATION_DEPTH:
#         run traverse_decomposition(handle) in a new process
#         return    # because the main process is only for dispatching the work to the child processes
#     else:
#         for child_handle in index.get_children(handle):
#             top_level_traversal_helper(child_handle, current_depth + 1)


# def top_level_traversal_before_parallelization(root_handle):
#     top_level_traversal_helper(root_handle, current_depth=0)


# # fetch root handle
# if name == "__main__":
#     root_handle = index.get_root()
#     top_level_traversal_before_parallelization(root_handle)
