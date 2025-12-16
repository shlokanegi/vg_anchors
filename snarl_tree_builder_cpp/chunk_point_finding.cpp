/**
 * Window Stats Generator and Chunk Point Finder for Variation Graph
 * 
 * This program analyzes a variation graph to:
 * 1. Compute sliding window statistics based on leaf snarl counts and haplotype informativeness
 * 2. Tag windows as "hap-dense" based on hap density threshold
 * 3. Find chunk points that ensure boundaries land in hap-dense windows
 * 
 * Chunk Point Finding Algorithm:
 * - Count leaf snarls until target is reached to get initial chunk boundary
 * - If boundary is in a hap-dense window, keep it
 * - If not, search forward and backward (respecting chain orientation) to find nearest hap-dense window
 * - Move boundary to that window
 * - Handle overlaps between consecutive chunks
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <thread>
#include <mutex>
#include <chrono>
#include <algorithm>
#include <memory>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <optional>

#include "bdsg/packed_graph.hpp"
#include "bdsg/snarl_distance_index.hpp"
#include "handlegraph/types.hpp"
#include "handlegraph/util.hpp"

using namespace std;
using namespace bdsg;
using namespace handlegraph;

// ============================================================================
// Data Structures
// ============================================================================

// Structure to hold the snarl tree map entry: (children_ids, leaf_count)
struct TreeMapEntry {
    vector<string> children_ids;
    int leaf_count;
};

// Structure to hold haplotype count info for a leaf snarl
struct HapCountInfo {
    string snarl_id;
    nid_t start_bound;
    nid_t end_bound;
    vector<int> step_counts;
    bool is_hap_informative;
};

// Structure to hold a leaf snarl with ordering info
struct LeafSnarlInfo {
    string snarl_id;
    nid_t start_node;
    nid_t end_node;
    bool is_hap_informative;
    size_t order_idx;
    string direct_child_id;  // ID of the direct child of top-level chain containing this leaf snarl
};

// Structure to hold window statistics
struct WindowStats {
    string chain_id;
    size_t window_idx;
    nid_t start_node;
    nid_t end_node;
    size_t leaf_snarl_count;
    size_t hap_informative_count;
    size_t window_length_bp;
    string path_name;
    double hap_density;         // hap_informative_count / window_length_bp * 1000 (per kb)
    bool is_hap_dense;          // true if hap_density >= threshold
    size_t start_leaf_idx;      // index of first leaf snarl in this window
    size_t end_leaf_idx;        // index of last leaf snarl in this window (exclusive)
};

// Structure to hold a chunk point
struct ChunkPoint {
    string chain_id;
    nid_t start_node;
    nid_t end_node;
    bool start_is_reverse;
    bool end_is_reverse;
    size_t leaf_snarl_count;
    size_t hap_informative_count;
    size_t chunk_length_bp;
    string path_name;
    size_t chunk_idx;           // sequential index within the chain
    size_t window_idx;          // window index where boundary was placed
    bool boundary_adjusted;     // true if boundary was moved from initial position
    vector<string> overlap_snarl_ids;
    size_t overlap_length_bp;
};

// ============================================================================
// Window Stats Builder and Chunk Finder Class
// ============================================================================

class ChunkPointFinder {
private:
    PackedGraph graph;
    SnarlDistanceIndex index;
    map<string, TreeMapEntry> snarl_tree_map;
    
    // Hap count data
    unordered_map<string, HapCountInfo> hap_counts;
    unordered_set<string> hap_informative_snarls;
    
    // Window parameters
    int window_size;
    int slide_size;
    int min_hap_support;
    double lower_percentile;      // lower percentile for threshold (e.g., 50 for median)
    double upper_percentile;      // upper percentile for threshold (e.g., 75 for 75th percentile)
    
    // Chunk parameters
    int target_leaf_snarls;
    int min_chunk_size;
    int overlap_snarl_count;
    
    // Results
    map<string, vector<WindowStats>> chain_windows;  // windows per chain
    map<string, vector<LeafSnarlInfo>> chain_leaf_snarls;  // leaf snarls per chain
    map<string, bool> chain_is_reverse;  // true if chain is traversed in reverse (decreasing node IDs)
    map<string, pair<double, double>> chr_thresholds;  // per-chromosome thresholds (min, max)
    vector<ChunkPoint> chunk_points;
    mutex results_mutex;
    
    // ========================================================================
    // Snarl ID Conversion Utilities
    // ========================================================================
    
    pair<nid_t, nid_t> parse_snarl_id(const string& snarl_id) {
        if (snarl_id.empty()) return {0, 0};
        
        size_t start_pos = (snarl_id[0] == 's' || snarl_id[0] == 'c') ? 1 : 0;
        size_t sep_pos = snarl_id.find('-', start_pos);
        if (sep_pos == string::npos) sep_pos = snarl_id.find('_', start_pos);
        if (sep_pos == string::npos || sep_pos <= start_pos) return {0, 0};
        
        try {
            nid_t node1 = stoull(snarl_id.substr(start_pos, sep_pos - start_pos));
            nid_t node2 = stoull(snarl_id.substr(sep_pos + 1));
            return {node1, node2};
        } catch (...) {
            return {0, 0};
        }
    }
    
    pair<string, string> json_to_tsv_keys(const string& json_snarl_id) {
        auto [node1, node2] = parse_snarl_id(json_snarl_id);
        if (node1 == 0 && node2 == 0) return {"", ""};
        return {to_string(node1) + "_" + to_string(node2), 
                to_string(node2) + "_" + to_string(node1)};
    }
    
    bool is_snarl_hap_informative(const string& json_snarl_id) {
        auto [key1, key2] = json_to_tsv_keys(json_snarl_id);
        return hap_informative_snarls.count(key1) > 0 || hap_informative_snarls.count(key2) > 0;
    }
    
    // ========================================================================
    // Graph/Index Utilities
    // ========================================================================
    
    pair<nid_t, nid_t> get_boundary_nodes(const net_handle_t& net) {
        net_handle_t left_bound = index.get_bound(net, false, false);
        net_handle_t right_bound = index.get_bound(net, true, false);
        nid_t left_id = graph.get_id(index.get_handle(left_bound, &graph));
        nid_t right_id = graph.get_id(index.get_handle(right_bound, &graph));
        return {left_id, right_id};
    }
    
    string get_id(const net_handle_t& net) {
        if (index.is_chain(net) || index.is_snarl(net)) {
            auto [left_id, right_id] = get_boundary_nodes(net);
            if (left_id < right_id) swap(left_id, right_id);
            string bound_str = to_string(left_id) + "-" + to_string(right_id);
            return (index.is_chain(net) ? "c" : "s") + bound_str;
        } else if (index.is_root(net)) {
            return "root";
        } else if (index.is_node(net)) {
            return to_string(graph.get_id(index.get_handle(net, &graph)));
        }
        throw runtime_error("Unknown net handle type");
    }
    
    string get_path_for_node(nid_t node_id) {
        string path_name = "unknown";
        handle_t node_handle = graph.get_handle(node_id);
        
        graph.for_each_step_on_handle(node_handle, [&](const step_handle_t& step) {
            path_handle_t path = graph.get_path_handle_of_step(step);
            string current_path_name = graph.get_path_name(path);
            
            if (current_path_name.find("CHM13") != string::npos || 
                current_path_name.find("chm13") != string::npos) {
                path_name = current_path_name;
                return false;
            }
            if (path_name == "unknown" && 
                (current_path_name.find("GRCh38") != string::npos || 
                 current_path_name.find("grch38") != string::npos)) {
                path_name = current_path_name;
            }
            if (path_name == "unknown") {
                path_name = current_path_name;
            }
            return true;
        });
        
        return path_name;
    }
    
    size_t calculate_distance(nid_t start_node, nid_t end_node) {
        if (start_node == end_node) {
            return graph.get_length(graph.get_handle(start_node));
        }
        
        size_t dist1 = index.minimum_distance(start_node, false, 0, end_node, false, 0, false, &graph);
        size_t dist2 = index.minimum_distance(start_node, true, 0, end_node, true, 0, false, &graph);
        size_t dist3 = index.minimum_distance(start_node, false, 0, end_node, true, 0, false, &graph);
        size_t dist4 = index.minimum_distance(start_node, true, 0, end_node, false, 0, false, &graph);
        
        size_t min_dist = SIZE_MAX;
        if (dist1 > 0 && dist1 < min_dist) min_dist = dist1;
        if (dist2 > 0 && dist2 < min_dist) min_dist = dist2;
        if (dist3 > 0 && dist3 < min_dist) min_dist = dist3;
        if (dist4 > 0 && dist4 < min_dist) min_dist = dist4;
        
        return (min_dist == SIZE_MAX) ? 0 : min_dist;
    }
    
    // ========================================================================
    // Leaf Snarl Collection
    // ========================================================================
    
    /**
     * Recursively collect leaf snarls, tracking which direct child of the top-level chain they belong to.
     * @param node_id Current node being traversed
     * @param leaf_snarls Output vector of leaf snarls
     * @param order_counter Counter for ordering leaf snarls
     * @param current_direct_child_id The direct child of the top-level chain that contains this subtree
     */
    void collect_leaf_snarls_recursive(const string& node_id, 
                                       vector<LeafSnarlInfo>& leaf_snarls,
                                       size_t& order_counter,
                                       const string& current_direct_child_id) {
        auto it = snarl_tree_map.find(node_id);
        if (it == snarl_tree_map.end()) return;
        
        const TreeMapEntry& entry = it->second;
        
        if (entry.children_ids.empty() && entry.leaf_count == 1) {
            auto [start, end] = parse_snarl_id(node_id);
            LeafSnarlInfo info;
            info.snarl_id = node_id;
            info.start_node = start;
            info.end_node = end;
            info.is_hap_informative = is_snarl_hap_informative(node_id);
            info.order_idx = order_counter++;
            info.direct_child_id = current_direct_child_id;
            leaf_snarls.push_back(info);
            return;
        }
        
        for (const string& child_id : entry.children_ids) {
            collect_leaf_snarls_recursive(child_id, leaf_snarls, order_counter, current_direct_child_id);
        }
    }
    
    /**
     * Collect all leaf snarls for a top-level chain, with each leaf snarl tagged
     * with the direct child of the chain that contains it.
     */
    vector<LeafSnarlInfo> collect_leaf_snarls_for_chain(const string& chain_id) {
        vector<LeafSnarlInfo> leaf_snarls;
        size_t order_counter = 0;
        
        // Get direct children of the chain
        auto it = snarl_tree_map.find(chain_id);
        if (it == snarl_tree_map.end()) return leaf_snarls;
        
        // For each direct child, collect its leaf snarls with the direct child ID tagged
        for (const string& direct_child_id : it->second.children_ids) {
            collect_leaf_snarls_recursive(direct_child_id, leaf_snarls, order_counter, direct_child_id);
        }
        
        return leaf_snarls;
    }
    
    // ========================================================================
    // Window Creation and Tagging
    // ========================================================================
    
    vector<WindowStats> create_windows_for_chain(const string& chain_id,
                                                  const vector<LeafSnarlInfo>& leaf_snarls) {
        vector<WindowStats> windows;
        if (leaf_snarls.empty()) return windows;
        
        auto [chain_start, chain_end] = parse_snarl_id(chain_id);
        string path_name = get_path_for_node(chain_start);
        
        size_t n = leaf_snarls.size();
        
        if (n <= (size_t)window_size) {
            WindowStats ws;
            ws.chain_id = chain_id;
            ws.window_idx = 0;
            ws.start_node = leaf_snarls.front().start_node;
            ws.end_node = leaf_snarls.back().end_node;
            ws.leaf_snarl_count = n;
            ws.hap_informative_count = 0;
            ws.start_leaf_idx = 0;
            ws.end_leaf_idx = n;
            for (const auto& ls : leaf_snarls) {
                if (ls.is_hap_informative) ws.hap_informative_count++;
            }
            ws.window_length_bp = calculate_distance(ws.start_node, ws.end_node);
            ws.path_name = path_name;
            ws.hap_density = (ws.window_length_bp > 0) ? 
                             (double)ws.hap_informative_count / ws.window_length_bp * 1000.0 : 0.0;
            // Window tagging will be done after all windows are created using chromosome-specific thresholds
            ws.is_hap_dense = false;  // Will be set later
            windows.push_back(ws);
            return windows;
        }
        
        size_t window_idx = 0;
        for (size_t start = 0; start < n; start += slide_size) {
            size_t end = min(start + window_size, n);
            
            if (end - start < (size_t)slide_size && start > 0) break;
            
            WindowStats ws;
            ws.chain_id = chain_id;
            ws.window_idx = window_idx++;
            ws.start_node = leaf_snarls[start].start_node;
            ws.end_node = leaf_snarls[end - 1].end_node;
            ws.leaf_snarl_count = end - start;
            ws.hap_informative_count = 0;
            ws.start_leaf_idx = start;
            ws.end_leaf_idx = end;
            
            for (size_t i = start; i < end; ++i) {
                if (leaf_snarls[i].is_hap_informative) ws.hap_informative_count++;
            }
            
            ws.window_length_bp = calculate_distance(ws.start_node, ws.end_node);
            ws.path_name = path_name;
            ws.hap_density = (ws.window_length_bp > 0) ? 
                             ((double)ws.hap_informative_count / ws.window_length_bp) * 1000.0 : 0.0;
            // Window tagging will be done after all windows are created using chromosome-specific thresholds
            ws.is_hap_dense = false;  // Will be set later
            windows.push_back(ws);
            
            if (end >= n) break;
        }
        
        return windows;
    }
    
    // ========================================================================
    // Find Nearest Hap-Dense Window
    // ========================================================================
    
    /**
     * Find the nearest hap-dense window to the given window index.
     * Considers chain orientation when determining "forward" and "backward".
     * Ensures selected windows don't start before chunk_start_leaf_idx.
     * 
     * @param windows              Vector of windows for this chain
     * @param current_idx          Current window index
     * @param is_reverse           True if chain is traversed in reverse
     * @param min_idx              Minimum window index (from previous chunk boundary)
     * @param chunk_start_leaf_idx Minimum leaf index that must be in the selected window
     * @return                     Index of nearest hap-dense window, or current_idx if none found
     */
    size_t find_nearest_hap_dense_window(const vector<WindowStats>& windows,
                                          size_t current_idx,
                                          bool is_reverse,
                                          size_t min_idx,
                                          size_t chunk_start_leaf_idx) {
        if (windows.empty() || current_idx >= windows.size()) {
            return current_idx;
        }
        
        // If current window is already hap-dense and doesn't start before chunk_start_leaf_idx, keep it
        if (windows[current_idx].is_hap_dense &&
            windows[current_idx].start_leaf_idx >= chunk_start_leaf_idx) {
            return current_idx;
        }
        
        // Search in both directions
        // "Forward" in chain traversal = higher indices if chain is forward, lower if reverse
        // "Backward" in chain traversal = lower indices if chain is forward, higher if reverse
        
        optional<size_t> forward_idx;
        optional<size_t> backward_idx;
        size_t forward_dist = SIZE_MAX;
        size_t backward_dist = SIZE_MAX;
        
        if (!is_reverse) {
            // Chain traversed forward: forward = higher indices, backward = lower indices
            // Search forward (higher indices)
            for (size_t i = current_idx + 1; i < windows.size(); ++i) {
                if (windows[i].is_hap_dense && 
                    windows[i].start_leaf_idx >= chunk_start_leaf_idx) {
                    forward_idx = i;
                    forward_dist = i - current_idx;
                    break;
                }
            }
            
            // Search backward (lower indices), but don't go below min_idx
            // AND ensure window doesn't start before chunk_start_leaf_idx
            for (size_t i = current_idx; i > min_idx; --i) {
                size_t check_idx = i - 1;
                if (windows[check_idx].is_hap_dense &&
                    windows[check_idx].start_leaf_idx >= chunk_start_leaf_idx) {
                    backward_idx = check_idx;
                    backward_dist = current_idx - check_idx;
                    break;
                }
            }
        } else {
            // Chain traversed in reverse: forward = lower indices, backward = higher indices
            // Search forward (lower indices), but don't go below min_idx
            // AND ensure window doesn't start before chunk_start_leaf_idx
            for (size_t i = current_idx; i > min_idx; --i) {
                size_t check_idx = i - 1;
                if (windows[check_idx].is_hap_dense &&
                    windows[check_idx].start_leaf_idx >= chunk_start_leaf_idx) {
                    forward_idx = check_idx;
                    forward_dist = current_idx - check_idx;
                    break;
                }
            }
            
            // Search backward (higher indices)
            for (size_t i = current_idx + 1; i < windows.size(); ++i) {
                if (windows[i].is_hap_dense &&
                    windows[i].start_leaf_idx >= chunk_start_leaf_idx) {
                    backward_idx = i;
                    backward_dist = i - current_idx;
                    break;
                }
            }
        }
        
        // Choose the nearest one, preferring forward on ties
        if (forward_idx.has_value() && backward_idx.has_value()) {
            if (forward_dist <= backward_dist) {
                return forward_idx.value();
            } else {
                return backward_idx.value();
            }
        } else if (forward_idx.has_value()) {
            return forward_idx.value();
        } else if (backward_idx.has_value()) {
            return backward_idx.value();
        }
        
        // No hap-dense window found, keep current
        return current_idx;
    }
    
    /**
     * Find which window contains a given leaf snarl index.
     */
    size_t find_window_for_leaf_idx(const vector<WindowStats>& windows, size_t leaf_idx) {
        for (size_t i = 0; i < windows.size(); ++i) {
            if (leaf_idx >= windows[i].start_leaf_idx && leaf_idx < windows[i].end_leaf_idx) {
                return i;
            }
        }
        // If not found in any window, return the last window
        return windows.empty() ? 0 : windows.size() - 1;
    }
    

    bool check_node_in_leaf_snarl(const LeafSnarlInfo& leaf_snarl, nid_t node_id) {
        if (node_id == leaf_snarl.start_node || node_id == leaf_snarl.end_node) {
            return true;
        }
        return false;
    }

    void print_d_to_b_leaf_snarl_count(const vector<LeafSnarlInfo>& leaf_snarls, nid_t d_node_id, nid_t b_node_id) {
        int cnt = 0;
        auto flag = false;
        for (size_t i = 0; i < leaf_snarls.size(); i++) {
            
            if (check_node_in_leaf_snarl(leaf_snarls[i], d_node_id)) {
                cout<<"found start (d node)" << d_node_id<<"in snarl: "<<leaf_snarls[i].snarl_id<<endl;
                flag = true;
                continue;
            }
            if (flag) {
                cnt++;
                if (check_node_in_leaf_snarl(leaf_snarls[i], b_node_id)) {
                    cout << "number of leaf snarls in leaf_snarl_list between d_node_id and b_node_id: " << cnt << endl;
                    exit(1);
                }
            }
        }
    }
    
    /**
     * Check if a node is a boundary of a given snarl/chain.
     */
    bool is_boundary_of(const string& snarl_id, nid_t node_id) {
        auto [start, end] = parse_snarl_id(snarl_id);
        return (node_id == start || node_id == end);
    }
    
    /**
     * Get boundary node for a chunk, ensuring it's a boundary of a direct child of the top-level chain.
     * Uses the pre-computed direct_child_id from leaf snarls for O(1) lookup.
     * 
     * @param chain_id The top-level chain ID
     * @param initial_node The initial boundary node from leaf snarl
     * @param is_at_chunk_start True if this is the start boundary of the chunk, false for end boundary
     * @param leaf_snarls Vector of all leaf snarls in the chain (with pre-computed direct_child_id)
     * @param leaf_idx Index of the leaf snarl that gave us initial_node
     * @return Boundary node that is a boundary of a direct child of the chain
     */
    nid_t get_boundary_node_on_direct_child(const string& chain_id, nid_t initial_node,
                                            bool is_at_chunk_start,
                                            const vector<LeafSnarlInfo>& leaf_snarls,
                                            size_t leaf_idx) {
        if (leaf_idx >= leaf_snarls.size()) {
            return initial_node;  // Fallback
        }
        
        // Get the direct child ID from the pre-computed value in the leaf snarl
        const string& direct_child_id = leaf_snarls[leaf_idx].direct_child_id;
        if (direct_child_id.empty()) {
            return initial_node;  // Fallback
        }
        
        // Check if the initial_node is already a boundary of this direct child
        if (is_boundary_of(direct_child_id, initial_node)) {
            return initial_node;
        }
        
        // Get boundaries of the direct child
        auto [child_node1, child_node2] = parse_snarl_id(direct_child_id);
        
        // Determine which boundary to use based on position in traversal
        // Use the leaf snarls to understand the traversal direction
        
        if (is_at_chunk_start) {
            // For start boundary: find the entry point of this direct child in traversal
            // Look at the first leaf snarl within this direct child
            // Find the first leaf snarl with this direct_child_id
            for (size_t i = 0; i < leaf_snarls.size(); ++i) {
                if (leaf_snarls[i].direct_child_id == direct_child_id) {
                    // First leaf snarl of this direct child - its start_node is the entry point
                    nid_t entry_node = leaf_snarls[i].start_node;
                    if (entry_node == child_node1 || entry_node == child_node2) {
                        return entry_node;
                    }
                    break;
                }
            }
            // Fallback: use the boundary closer to the leaf snarl start
            if (leaf_snarls[leaf_idx].start_node == child_node1 || 
                leaf_snarls[leaf_idx].start_node == child_node2) {
                return leaf_snarls[leaf_idx].start_node;
            }
            return child_node1;
        } else {
            // For end boundary: find the exit point of this direct child in traversal
            // Look at the last leaf snarl within this direct child
            // Find the last leaf snarl with this direct_child_id
            for (size_t i = leaf_snarls.size(); i > 0; --i) {
                if (leaf_snarls[i-1].direct_child_id == direct_child_id) {
                    // Last leaf snarl of this direct child - its end_node is the exit point
                    nid_t exit_node = leaf_snarls[i-1].end_node;
                    if (exit_node == child_node1 || exit_node == child_node2) {
                        return exit_node;
                    }
                    break;
                }
            }
            // Fallback: use the boundary closer to the leaf snarl end
            if (leaf_snarls[leaf_idx].end_node == child_node1 || 
                leaf_snarls[leaf_idx].end_node == child_node2) {
                return leaf_snarls[leaf_idx].end_node;
            }
            return child_node2;
        }
    }
    // ========================================================================
    // Chunk Point Generation
    // ========================================================================
    
    void generate_chunks_for_chain(const string& chain_id,
                                   const vector<LeafSnarlInfo>& leaf_snarls,
                                   const vector<WindowStats>& windows,
                                   bool is_reverse,
                                   vector<ChunkPoint>& local_chunks) {
        if (leaf_snarls.empty() || windows.empty()) return;
        
        auto [chain_start, chain_end] = parse_snarl_id(chain_id);
        string path_name = get_path_for_node(chain_start);
        
        size_t n = leaf_snarls.size();
        
        // DEBUG code
        // if (path_name == "CHM13#0#chr1") {
        //     for (size_t i = 1; i < n; i++) {
        //         if (leaf_snarls[i].start_node > leaf_snarls[i-1].end_node) {
        //             cout << "chain_id: " << chain_id << ": " << "leaf_snarls[i].start_node: " << leaf_snarls[i].start_node << " " << "leaf_snarls[i-1].end_node: " << leaf_snarls[i-1].end_node << " Error: leaf snarls are not sorted" << endl;
        //             exit(1);
        //         }
        //     }
        // }
        // nid_t d_node_id = stoull("c");
        // cout<<"d_node_id: " << d_node_id << endl;
        // nid_t b_node_id = stoull("85875154");
        // cout<<"b_node_id: " << b_node_id << endl;
        // print_d_to_b_leaf_snarl_count(leaf_snarls, d_node_id, b_node_id);
        // print_d_to_b_leaf_snarl_count(leaf_snarls, b_node_id, d_node_id);

        // Determine the fixed chain boundary nodes based on orientation
        // chain_id format: cLARGER-SMALLER, so chain_start > chain_end always
        // If is_reverse: we traverse from chain_start to chain_end (decreasing node IDs)
        // If !is_reverse: we traverse from chain_end to chain_start (increasing node IDs)
        nid_t fixed_chain_start = is_reverse ? chain_start : chain_end;
        nid_t fixed_chain_end = is_reverse ? chain_end : chain_start;
        
        // If chain is small, create a single chunk with fixed boundaries
        if ((int)n <= target_leaf_snarls) {
            ChunkPoint cp;
            cp.chain_id = chain_id;
            cp.start_node = fixed_chain_start;
            cp.end_node = fixed_chain_end;
            cp.start_is_reverse = is_reverse;
            cp.end_is_reverse = is_reverse;
            cp.leaf_snarl_count = n;
            cp.hap_informative_count = 0;
            for (const auto& ls : leaf_snarls) {
                if (ls.is_hap_informative) cp.hap_informative_count++;
            }
            cp.chunk_length_bp = calculate_distance(cp.start_node, cp.end_node);
            cp.path_name = path_name;
            cp.chunk_idx = 0;  // Single chunk
            cp.window_idx = 0;
            cp.boundary_adjusted = false;
            // First (and only) chunk: no overlap
            cp.overlap_snarl_ids.clear();
            cp.overlap_length_bp = 0;
            local_chunks.push_back(cp);
            return;
        }
        
        // Track overlap info to pass to the NEXT chunk
        vector<string> pending_overlap_snarl_ids;
        size_t pending_overlap_length_bp = 0;
        
        // Find chunk boundaries
        size_t chunk_start_leaf_idx = 0;
        size_t prev_boundary_window_idx = 0;
        bool is_first_chunk = true;
        size_t chunk_count = 0;
        
        while (chunk_start_leaf_idx < n) {
            // Initial boundary: count target_leaf_snarls from chunk start
            size_t initial_end_leaf_idx = min(chunk_start_leaf_idx + target_leaf_snarls, n);
            
            // Check if this will be the last chunk
            size_t remaining_after = n - initial_end_leaf_idx;
            bool will_be_last_chunk = (remaining_after == 0) || (remaining_after < (size_t)min_chunk_size);
            
            size_t actual_end_leaf_idx;
            size_t adjusted_window_idx = 0;
            bool boundary_adjusted = false;
            
            if (will_be_last_chunk) {
                // Last chunk: fix end to chain boundary, no adjustment
                actual_end_leaf_idx = n;
                adjusted_window_idx = windows.size() > 0 ? windows.size() - 1 : 0;
                boundary_adjusted = false;
            } else {
                // Not the last chunk: find hap-dense window for boundary
                size_t initial_window_idx = find_window_for_leaf_idx(windows, initial_end_leaf_idx - 1);
                
                // Find nearest hap-dense window
                // Ensure it doesn't start before chunk_start_leaf_idx to prevent chunks containing other chunks
                adjusted_window_idx = find_nearest_hap_dense_window(
                    windows, initial_window_idx, is_reverse, prev_boundary_window_idx, chunk_start_leaf_idx);
                
                // // temporarily disabling hap-informative boundary adjustment for simulating original chunk finder behaviour
                // adjusted_window_idx = initial_window_idx;

                boundary_adjusted = (adjusted_window_idx != initial_window_idx);
                
                // Determine the actual boundary leaf index from the adjusted window
                // find_nearest_hap_dense_window ensures the window doesn't start before chunk_start_leaf_idx
                if (adjusted_window_idx < windows.size()) {
                    actual_end_leaf_idx = windows[adjusted_window_idx].end_leaf_idx;
                } else {
                    actual_end_leaf_idx = n;
                }
                
                // Ensure we don't exceed the chain
                actual_end_leaf_idx = min(actual_end_leaf_idx, n);
                
                // CRITICAL: Ensure the boundary doesn't go backward from chunk start
                // This is a safety check in case window.end_leaf_idx somehow is before chunk_start_leaf_idx
                actual_end_leaf_idx = max(actual_end_leaf_idx, chunk_start_leaf_idx + 1);
                
                // Ensure minimum chunk size
                if (actual_end_leaf_idx < chunk_start_leaf_idx + (size_t)min_chunk_size && 
                    actual_end_leaf_idx < n) {
                    actual_end_leaf_idx = min(chunk_start_leaf_idx + min_chunk_size, n);
                }
                
                // Re-check if this becomes the last chunk
                remaining_after = n - actual_end_leaf_idx;
                if (remaining_after > 0 && remaining_after < (size_t)min_chunk_size) {
                    actual_end_leaf_idx = n;
                    will_be_last_chunk = true;
                }
            }
            
            // Create chunk
            ChunkPoint cp;
            cp.chain_id = chain_id;
            
            // Set start boundary
            if (is_first_chunk) {
                // First chunk: ALWAYS use exact chain start boundary
                // This ensures end-to-end coverage of the chain
                cp.start_node = fixed_chain_start;
            } else {
                // Intermediate chunk: get boundary from leaf snarl, adjusted to direct child boundary
                nid_t start_candidate = leaf_snarls[chunk_start_leaf_idx].start_node;
                cp.start_node = get_boundary_node_on_direct_child(
                    chain_id, start_candidate, true, leaf_snarls, chunk_start_leaf_idx);
            }
            
            // Set end boundary
            if (actual_end_leaf_idx >= n) {
                // Last chunk: ALWAYS use exact chain end boundary
                // This ensures end-to-end coverage of the chain
                cp.end_node = fixed_chain_end;
            } else {
                // Intermediate chunk: get boundary from leaf snarl, adjusted to direct child boundary
                nid_t end_candidate = leaf_snarls[actual_end_leaf_idx - 1].end_node;
                cp.end_node = get_boundary_node_on_direct_child(
                    chain_id, end_candidate, false, leaf_snarls, actual_end_leaf_idx - 1);
            }
            
            cp.start_is_reverse = is_reverse;
            cp.end_is_reverse = is_reverse;
            cp.leaf_snarl_count = actual_end_leaf_idx - chunk_start_leaf_idx;
            cp.hap_informative_count = 0;
            for (size_t i = chunk_start_leaf_idx; i < actual_end_leaf_idx; ++i) {
                if (leaf_snarls[i].is_hap_informative) cp.hap_informative_count++;
            }
            cp.chunk_length_bp = calculate_distance(cp.start_node, cp.end_node);
            cp.path_name = path_name;
            cp.chunk_idx = chunk_count;
            cp.window_idx = adjusted_window_idx;
            cp.boundary_adjusted = boundary_adjusted;
            
            // First chunk has NO overlap (nothing to overlap from)
            // Subsequent chunks receive the overlap from the previous chunk
            if (is_first_chunk) {
                cp.overlap_snarl_ids.clear();
                cp.overlap_length_bp = 0;
            } else {
                cp.overlap_snarl_ids = pending_overlap_snarl_ids;
                cp.overlap_length_bp = pending_overlap_length_bp;
            }
            
            local_chunks.push_back(cp);
            
            // Prepare overlap for the NEXT chunk (if this is not the last chunk)
            pending_overlap_snarl_ids.clear();
            pending_overlap_length_bp = 0;
            
            if (actual_end_leaf_idx < n && overlap_snarl_count > 0) {
                size_t overlap_start = (actual_end_leaf_idx > (size_t)overlap_snarl_count) ? 
                                       actual_end_leaf_idx - overlap_snarl_count : 0;
                // Don't include snarls before the current chunk start
                overlap_start = max(overlap_start, chunk_start_leaf_idx);
                
                for (size_t i = overlap_start; i < actual_end_leaf_idx; ++i) {
                    pending_overlap_snarl_ids.push_back(leaf_snarls[i].snarl_id);
                }
                if (!pending_overlap_snarl_ids.empty()) {
                    nid_t overlap_start_node = leaf_snarls[overlap_start].start_node;
                    nid_t overlap_end_node = leaf_snarls[actual_end_leaf_idx - 1].end_node;
                    pending_overlap_length_bp = calculate_distance(overlap_start_node, overlap_end_node);
                }
            }
            
            // Move to next chunk
            if (actual_end_leaf_idx >= n) break;
            
            // Next chunk starts from overlap region
            size_t next_start = (actual_end_leaf_idx > (size_t)overlap_snarl_count) ?
                               actual_end_leaf_idx - overlap_snarl_count : 0;
            
            // But don't go backward from current start
            next_start = max(next_start, chunk_start_leaf_idx + 1);
            
            chunk_start_leaf_idx = next_start;
            prev_boundary_window_idx = adjusted_window_idx;
            is_first_chunk = false;
            chunk_count++;
        }
    }
    
    // ========================================================================
    // Process Chain Worker
    // ========================================================================
    
    void process_chain(const string& chain_id) {
        // Collect leaf snarls
        vector<LeafSnarlInfo> leaf_snarls = collect_leaf_snarls_for_chain(chain_id);
        if (leaf_snarls.empty()) return;
        
        // Determine chain orientation
        nid_t first_node = leaf_snarls.front().start_node;
        nid_t last_node = leaf_snarls.back().end_node;
        bool is_reverse = (first_node > last_node);
        
        // Create windows
        vector<WindowStats> windows = create_windows_for_chain(chain_id, leaf_snarls);
        
        // Store results thread-safely
        {
            lock_guard<mutex> lock(results_mutex);
            chain_windows[chain_id] = windows;
            chain_leaf_snarls[chain_id] = leaf_snarls;
            chain_is_reverse[chain_id] = is_reverse;
        }
        
        // Chunk generation will be done after all windows are processed and thresholds are calculated
    }
    
    // ========================================================================
    // Chromosome-Specific Threshold Calculation
    // ========================================================================
    
    double percentile_value(const vector<double>& sorted_values, double percentile) {
        if (sorted_values.empty()) return 0.0;
        if (sorted_values.size() == 1) return sorted_values[0];
        
        // Clamp percentile to [0, 100]
        percentile = max(0.0, min(100.0, percentile));
        
        double position = (percentile / 100.0) * (sorted_values.size() - 1);
        size_t lower_idx = (size_t)floor(position);
        size_t upper_idx = (size_t)ceil(position);
        
        if (lower_idx == upper_idx) {
            return sorted_values[lower_idx];
        }
        
        double weight = position - lower_idx;
        return sorted_values[lower_idx] * (1.0 - weight) + sorted_values[upper_idx] * weight;
    }
    
    void calculate_chromosome_thresholds() {
        cerr << "\nCalculating chromosome-specific hap density thresholds..." << endl;
        cerr << "Using percentile range: [" << fixed << setprecision(1) << lower_percentile 
             << ", " << upper_percentile << "]" << endl;
        
        // Group windows by chromosome
        map<string, vector<double>> chr_densities;
        
        for (const auto& [chain_id, windows] : chain_windows) {
            for (const auto& ws : windows) {
                chr_densities[ws.path_name].push_back(ws.hap_density);
            }
        }
        
        // Calculate thresholds for each chromosome using specified percentiles
        for (auto& [chr, densities] : chr_densities) {
            if (densities.empty()) continue;
            
            sort(densities.begin(), densities.end());
            
            double lower_thresh = percentile_value(densities, lower_percentile);
            double upper_thresh = percentile_value(densities, upper_percentile);
            
            chr_thresholds[chr] = make_pair(lower_thresh, upper_thresh);
            
            cerr << "  " << chr << ": " << fixed << setprecision(1) << lower_percentile 
                 << "th=" << setprecision(3) << lower_thresh 
                 << ", " << setprecision(1) << upper_percentile 
                 << "th=" << setprecision(3) << upper_thresh << " per kb" << endl;
        }
    }
    
    void apply_chromosome_thresholds() {
        cerr << "Applying chromosome-specific thresholds to windows..." << endl;
        
        // Re-tag all windows with chromosome-specific thresholds
        for (auto& [chain_id, windows] : chain_windows) {
            for (auto& ws : windows) {
                auto it = chr_thresholds.find(ws.path_name);
                if (it != chr_thresholds.end()) {
                    double min_thresh = it->second.first;
                    double max_thresh = it->second.second;
                    ws.is_hap_dense = (ws.hap_density >= min_thresh) && 
                                      (ws.hap_density <= max_thresh);
                } else {
                    ws.is_hap_dense = false;
                }
            }
        }
    }
    
    void generate_chunks_from_windows(int num_threads = 0) {
        cerr << "Generating chunks using chromosome-specific thresholds..." << endl;
        
        chunk_points.clear();
        
        vector<string> chain_ids = get_top_level_chains();
        
        unsigned int actual_threads;
        if (num_threads > 0) {
            actual_threads = num_threads;
        } else {
            actual_threads = min((unsigned int)thread::hardware_concurrency(), 
                                 (unsigned int)chain_ids.size());
        }
        
        cerr << "Using " << actual_threads << " threads..." << endl;
        
        vector<thread> workers;
        size_t chains_per_thread = (chain_ids.size() + actual_threads - 1) / actual_threads;
        
        for (unsigned int t = 0; t < actual_threads; ++t) {
            workers.emplace_back([this, &chain_ids, t, chains_per_thread]() {
                size_t start = t * chains_per_thread;
                size_t end = min(start + chains_per_thread, chain_ids.size());
                
                for (size_t i = start; i < end; ++i) {
                    const string& chain_id = chain_ids[i];
                    
                    // Get stored windows and leaf snarls
                    vector<LeafSnarlInfo> leaf_snarls;
                    vector<WindowStats> windows;
                    bool is_reverse = false;
                    
                    {
                        lock_guard<mutex> lock(results_mutex);
                        auto it = chain_leaf_snarls.find(chain_id);
                        if (it != chain_leaf_snarls.end()) {
                            leaf_snarls = it->second;
                        }
                        auto win_it = chain_windows.find(chain_id);
                        if (win_it != chain_windows.end()) {
                            windows = win_it->second;
                        }
                        auto rev_it = chain_is_reverse.find(chain_id);
                        if (rev_it != chain_is_reverse.end()) {
                            is_reverse = rev_it->second;
                        }
                    }
                    
                    if (leaf_snarls.empty() || windows.empty()) continue;
                    
                    // Generate chunks
                    vector<ChunkPoint> local_chunks;
                    generate_chunks_for_chain(chain_id, leaf_snarls, windows, is_reverse, local_chunks);
                    
                    {
                        lock_guard<mutex> lock(results_mutex);
                        chunk_points.insert(chunk_points.end(), local_chunks.begin(), local_chunks.end());
                    }
                }
            });
        }
        
        for (auto& worker : workers) {
            worker.join();
        }
        
        auto end = chrono::high_resolution_clock::now();
        cerr << "Generated " << chunk_points.size() << " chunk points" << endl;
    }
    
public:
    ChunkPointFinder() 
        : window_size(1000), slide_size(500), min_hap_support(5),
          lower_percentile(50.0), upper_percentile(75.0),  // default: median to 75th percentile
          target_leaf_snarls(80000), min_chunk_size(35000), overlap_snarl_count(20) {}
    
    // Setters
    void set_window_size(int size) { window_size = size; }
    void set_slide_size(int size) { slide_size = size; }
    void set_min_hap_support(int support) { min_hap_support = support; }
    void set_lower_percentile(double pct) { lower_percentile = pct; }
    void set_upper_percentile(double pct) { upper_percentile = pct; }
    void set_target_leaf_snarls(int target) { target_leaf_snarls = target; }
    void set_min_chunk_size(int size) { min_chunk_size = size; }
    void set_overlap_snarl_count(int count) { overlap_snarl_count = count; }
    
    // ========================================================================
    // Data Loading
    // ========================================================================
    
    void load_graph_and_index(const string& graph_path, const string& index_path) {
        auto t0 = chrono::high_resolution_clock::now();
        
        ifstream graph_file(graph_path);
        if (!graph_file.is_open()) {
            throw runtime_error("Cannot open graph file: " + graph_path);
        }
        graph.deserialize(graph_file);
        graph_file.close();
        
        auto t1 = chrono::high_resolution_clock::now();
        cerr << "Loading graph: " << fixed << setprecision(2) 
             << chrono::duration<double>(t1 - t0).count() << " seconds" << endl;
        
        ifstream index_file(index_path);
        if (!index_file.is_open()) {
            throw runtime_error("Cannot open index file: " + index_path);
        }
        index.deserialize(index_file);
        index_file.close();
        
        auto t2 = chrono::high_resolution_clock::now();
        cerr << "Loading index: " << fixed << setprecision(2) 
             << chrono::duration<double>(t2 - t1).count() << " seconds" << endl;
    }
    
    void load_snarl_tree_json(const string& json_path) {
        auto start = chrono::high_resolution_clock::now();
        
        ifstream in(json_path);
        if (!in.is_open()) {
            throw runtime_error("Cannot open snarl tree JSON: " + json_path);
        }
        
        snarl_tree_map.clear();
        
        string line;
        string current_key;
        vector<string> children_ids;
        int leaf_count = 0;
        int state = 0;
        
        while (getline(in, line)) {
            if (line.empty() || line.find_first_not_of(" \t") == string::npos) continue;
            if (line.find("{") != string::npos && line.find("\"") == string::npos) continue;
            if (line.find("}") != string::npos && line.find("\"") == string::npos) continue;
            
            size_t key_start = line.find("\"");
            if (key_start != string::npos && state == 0) {
                size_t key_end = line.find("\"", key_start + 1);
                current_key = line.substr(key_start + 1, key_end - key_start - 1);
                state = 1;
                children_ids.clear();
                continue;
            }
            
            if (state == 1) {
                size_t child_start = line.find("\"");
                if (child_start != string::npos) {
                    size_t child_end = line.find("\"", child_start + 1);
                    string child_id = line.substr(child_start + 1, child_end - child_start - 1);
                    children_ids.push_back(child_id);
                } else if (line.find("]") != string::npos) {
                    state = 2;
                }
                continue;
            }
            
            if (state == 2) {
                size_t num_start = line.find_first_of("0123456789");
                if (num_start != string::npos) {
                    size_t num_end = line.find_first_not_of("0123456789", num_start);
                    string num_str = line.substr(num_start, num_end - num_start);
                    leaf_count = stoi(num_str);
                    snarl_tree_map[current_key] = {children_ids, leaf_count};
                    state = 0;
                }
            }
        }
        
        in.close();
        
        auto end = chrono::high_resolution_clock::now();
        cerr << "Loaded snarl tree: " << fixed << setprecision(2) 
             << chrono::duration<double>(end - start).count() << " seconds (" 
             << snarl_tree_map.size() << " entries)" << endl;
    }
    
    void load_hap_counts_tsv(const string& tsv_path) {
        auto start = chrono::high_resolution_clock::now();
        
        ifstream in(tsv_path);
        if (!in.is_open()) {
            throw runtime_error("Cannot open hap_counts TSV: " + tsv_path);
        }
        
        hap_counts.clear();
        hap_informative_snarls.clear();
        
        string line;
        bool header_skipped = false;
        size_t total_snarls = 0;
        size_t informative_snarls = 0;
        
        while (getline(in, line)) {
            if (line.empty()) continue;
            
            if (!header_skipped) {
                if (line.find("snarl_id") != string::npos) {
                    header_skipped = true;
                    continue;
                }
            }
            
            istringstream iss(line);
            string snarl_id, start_str, end_str, step_counts_str;
            
            if (!getline(iss, snarl_id, '\t')) continue;
            if (!getline(iss, start_str, '\t')) continue;
            if (!getline(iss, end_str, '\t')) continue;
            if (!getline(iss, step_counts_str, '\t')) continue;
            
            HapCountInfo info;
            info.snarl_id = snarl_id;
            
            try {
                info.start_bound = stoull(start_str);
                info.end_bound = stoull(end_str);
            } catch (...) {
                continue;
            }
            
            istringstream counts_stream(step_counts_str);
            string count_str;
            while (getline(counts_stream, count_str, ',')) {
                try {
                    info.step_counts.push_back(stoi(count_str));
                } catch (...) {}
            }
            
            int alleles_with_support = 0;
            for (int count : info.step_counts) {
                if (count >= min_hap_support) alleles_with_support++;
            }
            info.is_hap_informative = (alleles_with_support >= 2);
            
            hap_counts[snarl_id] = info;
            
            string reversed_key = to_string(info.end_bound) + "_" + to_string(info.start_bound);
            if (reversed_key != snarl_id) {
                hap_counts[reversed_key] = info;
            }
            
            if (info.is_hap_informative) {
                hap_informative_snarls.insert(snarl_id);
                hap_informative_snarls.insert(reversed_key);
                informative_snarls++;
            }
            
            total_snarls++;
        }
        
        in.close();
        
        auto end = chrono::high_resolution_clock::now();
        cerr << "Loaded hap_counts: " << fixed << setprecision(2) 
             << chrono::duration<double>(end - start).count() << " seconds (" 
             << total_snarls << " snarls, " << informative_snarls << " hap-informative)" << endl;
    }
    
    // ========================================================================
    // Processing
    // ========================================================================
    
    vector<string> get_top_level_chains() {
        vector<string> chains;
        auto root_it = snarl_tree_map.find("root");
        if (root_it != snarl_tree_map.end()) {
            for (const string& child_id : root_it->second.children_ids) {
                if (!child_id.empty() && child_id[0] == 'c') {
                    chains.push_back(child_id);
                }
            }
        }
        return chains;
    }
    
    void generate_chunks(int num_threads = 0) {
        auto start = chrono::high_resolution_clock::now();
        
        chunk_points.clear();
        chain_windows.clear();
        chain_leaf_snarls.clear();
        chain_is_reverse.clear();
        chr_thresholds.clear();
        
        vector<string> chain_ids = get_top_level_chains();
        cerr << "Processing " << chain_ids.size() << " top-level chains..." << endl;
        
        unsigned int actual_threads;
        if (num_threads > 0) {
            actual_threads = num_threads;
        } else {
            actual_threads = min((unsigned int)thread::hardware_concurrency(), 
                                (unsigned int)chain_ids.size());
        }
        if (actual_threads == 0) actual_threads = 1;
        
        cerr << "Using " << actual_threads << " threads..." << endl;
        
        vector<thread> threads;
        size_t chains_per_thread = (chain_ids.size() + actual_threads - 1) / actual_threads;
        
        for (size_t t = 0; t < actual_threads; ++t) {
            threads.emplace_back([this, t, &chain_ids, chains_per_thread]() {
                size_t start_idx = t * chains_per_thread;
                size_t end_idx = min(start_idx + chains_per_thread, chain_ids.size());
                
                for (size_t i = start_idx; i < end_idx; ++i) {
                    // leaf snarl collection and window creation inside this function
                    process_chain(chain_ids[i]);
                }
            });
        }
        
        for (auto& th : threads) {
            th.join();
        }
        
        // Calculate chromosome-specific hap-density thresholds, apply them, and generate chunks
        calculate_chromosome_thresholds();
        apply_chromosome_thresholds();
        generate_chunks_from_windows(num_threads);
        
        auto end = chrono::high_resolution_clock::now();
        cerr << "Generated " << chunk_points.size() << " chunk points in " 
             << fixed << setprecision(2) << chrono::duration<double>(end - start).count() 
             << " seconds" << endl;
    }
    
    // ========================================================================
    // Output
    // ========================================================================
    
    void write_window_stats(const string& output_path) {
        ofstream out(output_path);
        if (!out.is_open()) {
            throw runtime_error("Cannot open output file: " + output_path);
        }
        
        out << "chain_id\twindow_idx\tstart_node\tend_node\tleaf_snarl_count\t"
            << "hap_informative_count\twindow_length_bp\thap_density\tis_hap_dense\tpath_name\n";
        
        for (const auto& [chain_id, windows] : chain_windows) {
            for (const auto& ws : windows) {
                out << ws.chain_id << "\t"
                    << ws.window_idx << "\t"
                    << ws.start_node << "\t"
                    << ws.end_node << "\t"
                    << ws.leaf_snarl_count << "\t"
                    << ws.hap_informative_count << "\t"
                    << ws.window_length_bp << "\t"
                    << fixed << setprecision(4) << ws.hap_density << "\t"
                    << (ws.is_hap_dense ? "yes" : "no") << "\t"
                    << ws.path_name << "\n";
            }
        }
        
        out.close();
        cerr << "Window stats written to: " << output_path << endl;
    }
    
    /**
     * Extract chromosome number from path_name like "CHM13#0#chr20" -> 20
     * Returns 999 for non-standard chromosomes (e.g., chrX -> 23, chrY -> 24, chrM -> 25)
     */
    static int extract_chr_number(const string& path_name) {
        size_t chr_pos = path_name.find("chr");
        if (chr_pos == string::npos) {
            // Try "Chr" or "CHR"
            chr_pos = path_name.find("Chr");
            if (chr_pos == string::npos) {
                chr_pos = path_name.find("CHR");
            }
        }
        
        if (chr_pos == string::npos) {
            return 999;  // No chromosome found
        }
        
        string chr_part = path_name.substr(chr_pos + 3);
        
        // Handle special chromosomes
        if (chr_part[0] == 'X' || chr_part[0] == 'x') return 23;
        if (chr_part[0] == 'Y' || chr_part[0] == 'y') return 24;
        if (chr_part[0] == 'M' || chr_part[0] == 'm') return 25;
        
        // Extract numeric part
        try {
            size_t end_pos = chr_part.find_first_not_of("0123456789");
            if (end_pos == 0) return 999;
            return stoi(chr_part.substr(0, end_pos));
        } catch (...) {
            return 999;
        }
    }
    
    void normalize_reverse_chains() {
        // Group chunks by chain_id
        map<string, vector<ChunkPoint>> chunks_by_chain;
        for (const auto& cp : chunk_points) {
            chunks_by_chain[cp.chain_id].push_back(cp);
        }
        
        // Normalize chunks for reverse chains
        chunk_points.clear();
        for (auto& [chain_id, chunks] : chunks_by_chain) {
            // Sort chunks by chunk_idx first to ensure correct order
            sort(chunks.begin(), chunks.end(), 
                 [](const ChunkPoint& a, const ChunkPoint& b) {
                     return a.chunk_idx < b.chunk_idx;
                 });
            
            // Check if this chain is reverse
            bool is_reverse = chain_is_reverse.count(chain_id) && chain_is_reverse[chain_id];
            
            if (is_reverse) {
                // Store original overlaps before reversing
                vector<vector<string>> original_overlaps;
                vector<size_t> original_overlap_lengths;
                for (const auto& cp : chunks) {
                    original_overlaps.push_back(cp.overlap_snarl_ids);
                    original_overlap_lengths.push_back(cp.overlap_length_bp);
                }
                
                // Reverse the order of chunks
                reverse(chunks.begin(), chunks.end());
                
                // For each chunk, flip nodes and set forward orientation
                for (size_t i = 0; i < chunks.size(); ++i) {
                    auto& cp = chunks[i];
                    // Swap start_node and end_node
                    swap(cp.start_node, cp.end_node);
                    // Set both orientations to forward (+)
                    cp.start_is_reverse = false;
                    cp.end_is_reverse = false;
                }
                
                // Reassign overlaps after reversal
                // Original chunk i becomes new chunk (N-1-i)
                // The overlap on original chunk i should go to new chunk (N-1-i+1) if it exists
                // But we need to think about it differently:
                // - New chunk 0 (was last chunk): no overlap (first chunk)
                // - New chunk i (was chunk N-1-i): should get overlap from original chunk (N-i) if it exists
                //   which is the overlap that was on original chunk (N-i), which is now new chunk (i-1)
                size_t n = chunks.size();
                for (size_t i = 0; i < n; ++i) {
                    if (i == 0) {
                        // First chunk (was last chunk): no overlap
                        chunks[i].overlap_snarl_ids.clear();
                        chunks[i].overlap_length_bp = 0;
                    } else {
                        // New chunk i was original chunk (n-1-i)
                        // It should get the overlap that was on original chunk (n-i)
                        // which is now new chunk (i-1)
                        size_t orig_idx = n - i;
                        if (orig_idx < original_overlaps.size()) {
                            chunks[i].overlap_snarl_ids = original_overlaps[orig_idx];
                            chunks[i].overlap_length_bp = original_overlap_lengths[orig_idx];
                        } else {
                            chunks[i].overlap_snarl_ids.clear();
                            chunks[i].overlap_length_bp = 0;
                        }
                    }
                }
            }
            
            // Re-assign chunk_idx to maintain sequential ordering
            for (size_t i = 0; i < chunks.size(); ++i) {
                chunks[i].chunk_idx = i;
            }
            
            // Add normalized chunks back
            chunk_points.insert(chunk_points.end(), chunks.begin(), chunks.end());
        }
    }
    
    void write_chunk_points(const string& output_path) {
        // Normalize reverse chains to forward order before writing
        normalize_reverse_chains();
        
        ofstream out(output_path);
        if (!out.is_open()) {
            throw runtime_error("Cannot open output file: " + output_path);
        }
        
        out << "chain_id\tstart_node\tstart_orientation\tend_node\tend_orientation\t"
            << "leaf_snarl_count\thap_informative_count\tchunk_length_bp\tpath_name\t"
            << "window_idx\tboundary_adjusted\toverlap_length_bp\toverlap_snarl_ids\n";
        
        // Sort by chromosome order, then by chunk order within each chain
        sort(chunk_points.begin(), chunk_points.end(), 
             [](const ChunkPoint& a, const ChunkPoint& b) {
                 int chr_a = extract_chr_number(a.path_name);
                 int chr_b = extract_chr_number(b.path_name);
                 if (chr_a != chr_b) return chr_a < chr_b;
                 if (a.chain_id != b.chain_id) return a.chain_id < b.chain_id;
                 return a.chunk_idx < b.chunk_idx;
             });
        
        for (const auto& cp : chunk_points) {
            out << cp.chain_id << "\t"
                << cp.start_node << "\t"
                << (cp.start_is_reverse ? "-" : "+") << "\t"
                << cp.end_node << "\t"
                << (cp.end_is_reverse ? "-" : "+") << "\t"
                << cp.leaf_snarl_count << "\t"
                << cp.hap_informative_count << "\t"
                << cp.chunk_length_bp << "\t"
                << cp.path_name << "\t"
                << cp.window_idx << "\t"
                << (cp.boundary_adjusted ? "yes" : "no") << "\t"
                << cp.overlap_length_bp << "\t";
            
            if (cp.overlap_snarl_ids.empty()) {
                out << "NA";
            } else {
                for (size_t i = 0; i < cp.overlap_snarl_ids.size(); ++i) {
                    out << cp.overlap_snarl_ids[i];
                    if (i < cp.overlap_snarl_ids.size() - 1) out << ",";
                }
            }
            out << "\n";
        }
        
        out.close();
        cerr << "Chunk points written to: " << output_path << endl;
    }
    
    void print_summary() {
        // Window statistics
        size_t total_windows = 0;
        size_t hap_dense_windows = 0;
        size_t total_leaf_snarls_in_chains = 0;
        size_t chains_with_windows = 0;
        
        // Per-chromosome statistics
        map<string, size_t> chr_windows;
        map<string, size_t> chr_hap_dense_windows;
        map<string, size_t> chr_leaf_snarls;
        map<string, size_t> chr_hap_informative;
        
        for (const auto& [chain_id, windows] : chain_windows) {
            total_windows += windows.size();
            if (!windows.empty()) chains_with_windows++;
            for (const auto& ws : windows) {
                if (ws.is_hap_dense) hap_dense_windows++;
                chr_windows[ws.path_name]++;
                if (ws.is_hap_dense) chr_hap_dense_windows[ws.path_name]++;
            }
        }
        
        for (const auto& [chain_id, leaf_snarls] : chain_leaf_snarls) {
            total_leaf_snarls_in_chains += leaf_snarls.size();
            if (!leaf_snarls.empty()) {
                string path_name = get_path_for_node(parse_snarl_id(chain_id).first);
                chr_leaf_snarls[path_name] += leaf_snarls.size();
                for (const auto& ls : leaf_snarls) {
                    if (ls.is_hap_informative) chr_hap_informative[path_name]++;
                }
            }
        }
        
        cerr << "\n=== Window Summary ===" << endl;
        cerr << "Chains processed: " << chain_windows.size() << endl;
        cerr << "Chains with windows: " << chains_with_windows << endl;
        cerr << "Total leaf snarls in chains: " << total_leaf_snarls_in_chains << endl;
        cerr << "Total windows: " << total_windows << endl;
        cerr << "Hap-dense windows: " << hap_dense_windows << " (" 
             << fixed << setprecision(1) << (100.0 * hap_dense_windows / max(total_windows, (size_t)1))
             << "%)" << endl;
        
        // Per-chromosome statistics
        cerr << "\n=== Per-Chromosome Statistics ===" << endl;
        cerr << left << setw(15) << "Chromosome" 
             << right << setw(12) << "Windows" 
             << setw(18) << "Hap-Dense Win" 
             << setw(15) << "Hap-Dense %"
             << setw(18) << "Leaf Snarls"
             << setw(20) << "Hap-Informative"
             << setw(12) << "Threshold" << endl;
        cerr << string(110, '-') << endl;
        
        // Sort chromosomes
        vector<string> sorted_chrs;
        for (const auto& [chr, _] : chr_windows) {
            sorted_chrs.push_back(chr);
        }
        sort(sorted_chrs.begin(), sorted_chrs.end(), [this](const string& a, const string& b) {
            return extract_chr_number(a) < extract_chr_number(b);
        });
        
        for (const string& chr : sorted_chrs) {
            size_t win_count = chr_windows[chr];
            size_t hap_dense_count = chr_hap_dense_windows[chr];
            size_t leaf_count = chr_leaf_snarls[chr];
            size_t hap_info_count = chr_hap_informative[chr];
            double hap_dense_pct = (win_count > 0) ? 100.0 * hap_dense_count / win_count : 0.0;
            
            string threshold_str = "N/A";
            auto it = chr_thresholds.find(chr);
            if (it != chr_thresholds.end()) {
                threshold_str = "[" + to_string(it->second.first) + ", " + 
                               to_string(it->second.second) + "]";
            }
            
            cerr << left << setw(15) << chr
                 << right << setw(12) << win_count
                 << setw(18) << hap_dense_count
                 << setw(15) << fixed << setprecision(1) << hap_dense_pct
                 << setw(18) << leaf_count
                 << setw(20) << hap_info_count
                 << setw(12) << threshold_str << endl;
        }
        
        // Chunk statistics
        size_t adjusted_boundaries = 0;
        size_t total_leaf_snarls = 0;
        size_t total_hap_informative = 0;
        for (const auto& cp : chunk_points) {
            if (cp.boundary_adjusted) adjusted_boundaries++;
            total_leaf_snarls += cp.leaf_snarl_count;
            total_hap_informative += cp.hap_informative_count;
        }
        
        cerr << "\n=== Chunk Summary ===" << endl;
        cerr << "Total chunks: " << chunk_points.size() << endl;
        cerr << "Adjusted boundaries: " << adjusted_boundaries << " (" 
             << fixed << setprecision(1) << (100.0 * adjusted_boundaries / max(chunk_points.size(), (size_t)1))
             << "%)" << endl;
        cerr << "Total leaf snarls: " << total_leaf_snarls << endl;
        cerr << "Total hap-informative: " << total_hap_informative << endl;
        if (chunk_points.size() > 0) {
            cerr << "Average chunk size: " << (total_leaf_snarls / chunk_points.size()) << " leaf snarls" << endl;
        }
    }
};

// ============================================================================
// Main Function
// ============================================================================

int main(int argc, char* argv[]) {
    string graph_path, index_path;
    string snarl_tree_json_path;
    string hap_counts_tsv_path;
    string window_stats_output = "";
    string chunk_points_output = "";
    int num_threads = 0;
    int window_size = 1000;
    int slide_size = 500;
    int min_hap_support = 5;
    double lower_percentile = 50.0;  // lower percentile (e.g., 50 for median)
    double upper_percentile = 75.0;  // upper percentile (e.g., 75 for 75th percentile)
    int target_leaf_snarls = 80000;
    int min_chunk_size = 35000;
    int overlap_snarl_count = 20;
    
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if ((arg == "-g" || arg == "--graph") && i + 1 < argc) {
            graph_path = argv[++i];
        } else if ((arg == "-i" || arg == "--index") && i + 1 < argc) {
            index_path = argv[++i];
        } else if ((arg == "-j" || arg == "--snarl-tree-json") && i + 1 < argc) {
            snarl_tree_json_path = argv[++i];
        } else if ((arg == "-H" || arg == "--hap-counts") && i + 1 < argc) {
            hap_counts_tsv_path = argv[++i];
        } else if ((arg == "-W" || arg == "--window-stats-output") && i + 1 < argc) {
            window_stats_output = argv[++i];
        } else if ((arg == "-c" || arg == "--chunk-points-output") && i + 1 < argc) {
            chunk_points_output = argv[++i];
        } else if ((arg == "-n" || arg == "--num-threads") && i + 1 < argc) {
            num_threads = stoi(argv[++i]);
        } else if ((arg == "-w" || arg == "--window-size") && i + 1 < argc) {
            window_size = stoi(argv[++i]);
        } else if ((arg == "-s" || arg == "--slide-size") && i + 1 < argc) {
            slide_size = stoi(argv[++i]);
        } else if ((arg == "-m" || arg == "--min-hap-support") && i + 1 < argc) {
            min_hap_support = stoi(argv[++i]);
        } else if ((arg == "-d" || arg == "--lower-percentile") && i + 1 < argc) {
            lower_percentile = stod(argv[++i]);
        } else if ((arg == "-D" || arg == "--upper-percentile") && i + 1 < argc) {
            upper_percentile = stod(argv[++i]);
        } else if ((arg == "-t" || arg == "--target-leaf-snarls") && i + 1 < argc) {
            target_leaf_snarls = stoi(argv[++i]);
        } else if ((arg == "-M" || arg == "--min-chunk-size") && i + 1 < argc) {
            min_chunk_size = stoi(argv[++i]);
        } else if ((arg == "-v" || arg == "--overlap-snarl-count") && i + 1 < argc) {
            overlap_snarl_count = stoi(argv[++i]);
        } else if (arg == "-h" || arg == "--help") {
            cout << "Usage: " << argv[0] << " [OPTIONS]\n";
            cout << "\nChunk Point Finder with Hap-Dense Window Constraints\n";
            cout << "====================================================\n";
            cout << "\nRequired Arguments:\n";
            cout << "  -g, --graph                  Path to variation graph (.pg/.vg)\n";
            cout << "  -i, --index                  Path to snarl distance index (.dist)\n";
            cout << "  -j, --snarl-tree-json        Path to pre-built snarl tree JSON\n";
            cout << "  -H, --hap-counts             Path to hap_counts TSV file\n";
            cout << "\nOutput Arguments (at least one required):\n";
            cout << "  -W, --window-stats-output    Output window stats TSV\n";
            cout << "  -c, --chunk-points-output    Output chunk points TSV\n";
            cout << "\nWindow Parameters:\n";
            cout << "  -w, --window-size            Leaf snarls per window (default: 1000)\n";
            cout << "  -s, --slide-size             Window slide size (default: 500)\n";
            cout << "  -m, --min-hap-support        Min haplotype support (default: 5)\n";
            cout << "  -d, --lower-percentile       Lower percentile for threshold (default: 50.0 = median)\n";
            cout << "  -D, --upper-percentile       Upper percentile for threshold (default: 75.0 = 75th percentile)\n";
            cout << "                               Thresholds are calculated per chromosome from hap density distribution\n";
            cout << "\nChunk Parameters:\n";
            cout << "  -t, --target-leaf-snarls     Target leaf snarls per chunk (default: 80000)\n";
            cout << "  -M, --min-chunk-size         Minimum chunk size (default: 35000)\n";
            cout << "  -v, --overlap-snarl-count    Min overlap snarls between chunks (default: 20)\n";
            cout << "\nOther:\n";
            cout << "  -n, --num-threads            Number of threads (default: auto)\n";
            cout << "  -h, --help                   Show this help\n";
            cout << "\nAlgorithm:\n";
            cout << "  1. Create sliding windows across each chain\n";
            cout << "  2. Tag windows as 'hap-dense' if min_threshold <= hap_density <= max_threshold\n";
            cout << "  3. Find initial chunk boundaries by counting target leaf snarls\n";
            cout << "  4. If boundary is not in hap-dense window, move to nearest one\n";
            cout << "     (respecting chain orientation for forward/backward search)\n";
            cout << "  5. Maintain overlap between consecutive chunks\n";
            cout << "\nExample:\n";
            cout << "  " << argv[0] << " \\\n";
            cout << "    -g graph.pg -i index.dist \\\n";
            cout << "    -j snarl_tree_map.json -H hap_counts.tsv \\\n";
            cout << "    -W window_stats.tsv -c chunk_points.tsv \\\n";
            cout << "    -d 1.0 -t 80000 -v 20 -n 32\n";
            return 0;
        }
    }
    
    // Validate arguments
    if (graph_path.empty() || index_path.empty() || 
        snarl_tree_json_path.empty() || hap_counts_tsv_path.empty()) {
        cerr << "Error: Missing required arguments.\n";
        cerr << "Required: -g GRAPH -i INDEX -j JSON -H HAP_COUNTS\n";
        cerr << "Use -h or --help for usage information.\n";
        return 1;
    }
    
    if (window_stats_output.empty() && chunk_points_output.empty()) {
        cerr << "Error: At least one output file required (-W or -c).\n";
        cerr << "Use -h or --help for usage information.\n";
        return 1;
    }
    
    try {
        ChunkPointFinder finder;
        finder.set_window_size(window_size);
        finder.set_slide_size(slide_size);
        finder.set_min_hap_support(min_hap_support);
        finder.set_lower_percentile(lower_percentile);
        finder.set_upper_percentile(upper_percentile);
        finder.set_target_leaf_snarls(target_leaf_snarls);
        finder.set_min_chunk_size(min_chunk_size);
        finder.set_overlap_snarl_count(overlap_snarl_count);
        
        // Load data
        finder.load_graph_and_index(graph_path, index_path);
        finder.load_snarl_tree_json(snarl_tree_json_path);
        finder.load_hap_counts_tsv(hap_counts_tsv_path);
        
        // Generate chunks (this also creates windows)
        finder.generate_chunks(num_threads);
        
        // Write outputs
        if (!window_stats_output.empty()) {
            finder.write_window_stats(window_stats_output);
        }
        if (!chunk_points_output.empty()) {
            finder.write_chunk_points(chunk_points_output);
        }
        
        finder.print_summary();
        
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}
