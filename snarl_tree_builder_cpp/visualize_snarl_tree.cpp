/**
 * Fast C++ implementation for visualizing snarl trees.
 * 
 * This program reads a snarl tree JSON and generates an HTML visualization
 * using multithreading for efficient processing of large trees.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <thread>
#include <mutex>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <memory>
#include <atomic>

// Simple JSON parser for our specific format
// We only need to parse: {"key": [[children...], count], ...}
class SnarlTreeData {
public:
    std::vector<std::string> children;
    int leaf_count;
    
    SnarlTreeData() : leaf_count(0) {}
    SnarlTreeData(const std::vector<std::string>& c, int lc) : children(c), leaf_count(lc) {}
};

class SnarlTreeMap {
private:
    std::map<std::string, SnarlTreeData> data;
    
public:
    void load_from_json(const std::string& filepath) {
        auto start = std::chrono::high_resolution_clock::now();
        std::cerr << "Loading JSON file..." << std::endl;
        
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filepath);
        }
        
        // Read entire file
        std::stringstream buffer;
        buffer << file.rdbuf();
        std::string content = buffer.str();
        file.close();
        
        // Parse JSON manually for speed
        parse_json(content);
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cerr << "Loaded " << data.size() << " entries in " 
                  << std::fixed << std::setprecision(2) << elapsed.count() 
                  << " seconds" << std::endl;
    }
    
    const SnarlTreeData* get(const std::string& key) const {
        auto it = data.find(key);
        if (it != data.end()) {
            return &it->second;
        }
        return nullptr;
    }
    
    size_t size() const { return data.size(); }
    
private:
    void parse_json(const std::string& content) {
        size_t pos = 1; // Skip opening '{'
        
        while (pos < content.length()) {
            // Skip whitespace
            while (pos < content.length() && std::isspace(content[pos])) pos++;
            
            if (content[pos] == '}') break;
            
            // Parse key
            if (content[pos] != '"') break;
            pos++; // Skip opening quote
            
            size_t key_start = pos;
            while (pos < content.length() && content[pos] != '"') pos++;
            std::string key = content.substr(key_start, pos - key_start);
            pos++; // Skip closing quote
            
            // Skip to '['
            while (pos < content.length() && content[pos] != '[') pos++;
            pos++; // Skip '['
            
            // Parse first array (children)
            while (pos < content.length() && content[pos] != '[') pos++;
            pos++; // Skip opening '['
            
            std::vector<std::string> children;
            while (pos < content.length() && content[pos] != ']') {
                while (pos < content.length() && std::isspace(content[pos])) pos++;
                if (content[pos] == '"') {
                    pos++; // Skip opening quote
                    size_t child_start = pos;
                    while (pos < content.length() && content[pos] != '"') pos++;
                    children.push_back(content.substr(child_start, pos - child_start));
                    pos++; // Skip closing quote
                }
                if (content[pos] == ',') pos++;
            }
            pos++; // Skip closing ']'
            
            // Skip to number
            while (pos < content.length() && (std::isspace(content[pos]) || content[pos] == ',')) pos++;
            
            // Parse leaf count
            size_t num_start = pos;
            while (pos < content.length() && std::isdigit(content[pos])) pos++;
            int leaf_count = std::stoi(content.substr(num_start, pos - num_start));
            
            data[key] = SnarlTreeData(children, leaf_count);
            
            // Skip to next entry
            while (pos < content.length() && content[pos] != ',' && content[pos] != '}') pos++;
            if (content[pos] == ',') pos++;
        }
    }
};

// Tree node for visualization
struct TreeNode {
    std::string name;
    std::string original_name;
    int num_leaf_snarls;
    std::vector<std::shared_ptr<TreeNode>> children;
    
    TreeNode(const std::string& n, const std::string& on, int nls)
        : name(n), original_name(on), num_leaf_snarls(nls) {}
};

class TreeBuilder {
private:
    const SnarlTreeMap& tree_map;
    std::atomic<int> chain_counter{1};
    std::atomic<int> snarl_counter{1};
    std::mutex output_mutex;
    
    std::string format_name(const std::string& name) {
        if (name == "root") return "Root";
        
        char type = name[0];
        std::string rest = name.substr(1);
        
        if (type == 'c') {
            return "Chain " + rest;
        } else if (type == 's') {
            return "Snarl " + rest;
        }
        return "Node " + name;
    }
    
    std::shared_ptr<TreeNode> build_tree_recursive(const std::string& node_id, int& local_chain, int& local_snarl) {
        std::string short_name;
        
        if (node_id == "root") {
            short_name = "Root";
        } else if (node_id[0] == 'c') {
            short_name = "Chain " + std::to_string(local_chain++);
        } else if (node_id[0] == 's') {
            short_name = "Snarl " + std::to_string(local_snarl++);
        }
        
        const SnarlTreeData* data = tree_map.get(node_id);
        
        if (!data) {
            // Leaf node not in map
            return std::make_shared<TreeNode>(short_name, format_name(node_id), 1);
        }
        
        auto node = std::make_shared<TreeNode>(short_name, format_name(node_id), data->leaf_count);
        
        for (const auto& child_id : data->children) {
            node->children.push_back(build_tree_recursive(child_id, local_chain, local_snarl));
        }
        
        return node;
    }
    
    void write_json_node(std::ostream& out, const TreeNode& node, int indent = 0) {
        std::string ind(indent * 4, ' ');
        
        out << ind << "{\n";
        out << ind << "    \"name\": \"" << json_escape(node.name) << "\",\n";
        out << ind << "    \"original_name\": \"" << json_escape(node.original_name) << "\",\n";
        out << ind << "    \"num_leaf_snarls\": " << node.num_leaf_snarls;
        
        if (!node.children.empty()) {
            out << ",\n" << ind << "    \"children\": [\n";
            for (size_t i = 0; i < node.children.size(); ++i) {
                write_json_node(out, *node.children[i], indent + 2);
                if (i < node.children.size() - 1) {
                    out << ",\n";
                } else {
                    out << "\n";
                }
            }
            out << ind << "    ]\n";
        } else {
            out << "\n";
        }
        
        out << ind << "}";
    }
    
    std::string json_escape(const std::string& s) {
        std::ostringstream o;
        for (char c : s) {
            switch (c) {
                case '"': o << "\\\""; break;
                case '\\': o << "\\\\"; break;
                case '\b': o << "\\b"; break;
                case '\f': o << "\\f"; break;
                case '\n': o << "\\n"; break;
                case '\r': o << "\\r"; break;
                case '\t': o << "\\t"; break;
                default:
                    if (c >= 0 && c <= 0x1f) {
                        o << "\\u" << std::hex << std::setw(4) << std::setfill('0') << (int)c;
                    } else {
                        o << c;
                    }
            }
        }
        return o.str();
    }
    
public:
    TreeBuilder(const SnarlTreeMap& tm) : tree_map(tm) {}
    
    std::shared_ptr<TreeNode> build_tree(const std::string& root_id = "root") {
        auto start = std::chrono::high_resolution_clock::now();
        std::cerr << "Building tree structure..." << std::endl;
        
        int chain_cnt = chain_counter.load();
        int snarl_cnt = snarl_counter.load();
        auto tree = build_tree_recursive(root_id, chain_cnt, snarl_cnt);
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cerr << "Tree built in " << std::fixed << std::setprecision(2) 
                  << elapsed.count() << " seconds" << std::endl;
        
        return tree;
    }
    
    void write_tree_json(std::ostream& out, const TreeNode& root) {
        write_json_node(out, root, 0);
    }
};

void write_html_template(const std::string& output_path, const std::string& tree_json,
                        int initial_depth, int filter_depth, int min_leaf_snarls) {
    std::cerr << "Writing HTML file..." << std::endl;
    
    std::ofstream out(output_path);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + output_path);
    }
    
    out << R"HTML(<!DOCTYPE html>
<html>
<head>
    <title>Snarl Tree Visualization</title>
    <style>
        .node circle {
            fill: #fff;
            stroke: steelblue;
            stroke-width: 1.5px;
        }
        .node--internal circle {
            cursor: pointer;
        }
        .node--leaf circle {
            fill: #fff;
        }
        .node text {
            font: 10px sans-serif;
            cursor: pointer;
        }
        .link {
            fill: none;
            stroke: #ccc;
            stroke-width: 1.5px;
        }
        #controls {
            position: absolute;
            top: 10px;
            left: 10px;
            font-family: sans-serif;
            background: rgba(240,240,240,0.8);
            padding: 10px;
            border: 1px solid #ccc;
            border-radius: 5px;
        }
    </style>
</head>
<body>
    <div id="controls">
        <div>
            <label for="depth-input">Display Depth:</label>
            <input type="number" id="depth-input" value=")HTML";
    out << initial_depth;
    out << R"HTML(" min="0" style="width: 50px;">
        </div>
        <div style="margin-top: 5px;">
            <label for="filter-depth-input">Filter Depth:</label>
            <input type="number" id="filter-depth-input" value=")HTML";
    out << filter_depth;
    out << R"HTML(" min="0" style="width: 50px;">
        </div>
        <div style="margin-top: 5px;">
            <label for="filter-input">Min Leaf Snarls:</label>
            <input type="number" id="filter-input" value=")HTML";
    out << min_leaf_snarls;
    out << R"HTML(" min="0" style="width: 50px;">
        </div>
        <button id="update-button" style="margin-top: 5px;">Update</button>
    </div>
    <div id="details" style="position: fixed; top: 10px; right: 10px; padding: 10px; background: rgba(240,240,240,0.8); border: 1px solid #ccc; border-radius: 5px; font-family: sans-serif;">
        Click a node to see details.
    </div>
    <script src="https://d3js.org/d3.v5.min.js"></script>
    <script>
        const fullTreeData = )HTML";
    
    out << tree_json;
    
    out << R"HTML(;
        let root;
        const margin = {top: 20, right: 120, bottom: 20, left: 120};
        const svg = d3.select("body").append("svg").append("g");
        let i = 0;
        const duration = 750;

        function initialize(treeData) {
            if (!treeData) {
                svg.selectAll("*").remove();
                return;
            }
            root = d3.hierarchy(treeData, d => d.children);
            root.x0 = 0;
            root.y0 = 0;
            const depth = parseInt(d3.select("#depth-input").property("value"));
            collapseAll(root);
            expandToDepth(root, depth);
            update(root);
        }

        d3.select("#update-button").on("click", () => {
            const minSnarls = parseInt(d3.select("#filter-input").property("value"));
            const filterDepth = parseInt(d3.select("#filter-depth-input").property("value"));
            if (!isNaN(minSnarls) && !isNaN(filterDepth)) {
                const treeDataCopy = JSON.parse(JSON.stringify(fullTreeData));
                const filteredData = filterData(treeDataCopy, minSnarls, filterDepth);
                initialize(filteredData);
            }
        });

        function filterData(node, minSnarls, filterDepth, currentDepth = 0) {
            if (!node) return null;
            let visibleChildren = [];
            if (node.children) {
                visibleChildren = node.children
                    .map(child => filterData(child, minSnarls, filterDepth, currentDepth + 1))
                    .filter(child => child !== null);
            }
            const hasVisibleChildren = visibleChildren.length > 0;
            if (currentDepth === filterDepth && node.num_leaf_snarls < minSnarls) return null;
            if (currentDepth < filterDepth && !hasVisibleChildren) return null;
            const newNode = { ...node };
            if (hasVisibleChildren) {
                newNode.children = visibleChildren;
            } else {
                delete newNode.children;
            }
            return newNode;
        }

        function collapseAll(d) {
            if (d.children) {
                d._children = d.children;
                d._children.forEach(collapseAll);
                d.children = null;
            }
        }

        function expandToDepth(d, depth) {
            if (depth > 0) {
                if (d._children) {
                    d.children = d._children;
                    d._children = null;
                }
                if (d.children) {
                    d.children.forEach(child => expandToDepth(child, depth - 1));
                }
            }
        }

        d3.select("#update-button").dispatch('click');

        function update(source) {
            const treeLayout = d3.tree().nodeSize([40, 180]);
            const treeData = treeLayout(root);
            const nodes = treeData.descendants();
            const links = treeData.descendants().slice(1);
            nodes.forEach(d => { d.y = d.depth * 250 });

            const node = svg.selectAll('g.node').data(nodes, d => d.id || (d.id = ++i));
            const nodeEnter = node.enter().append('g')
                .attr('class', 'node')
                .attr("transform", d => "translate(" + source.y0 + "," + source.x0 + ")")
                .on('click', click);

            nodeEnter.append('circle')
                .attr('class', 'node')
                .attr('r', 1e-6)
                .style("fill", d => d._children ? "lightsteelblue" : "#fff");

            nodeEnter.append('text')
                .attr("dy", ".35em")
                .attr("x", d => d.children || d._children ? -13 : 13)
                .attr("text-anchor", d => d.children || d._children ? "end" : "start")
                .text(d => `${d.data.name} (${d.data.num_leaf_snarls})`);

            const nodeUpdate = nodeEnter.merge(node);
            nodeUpdate.transition().duration(duration)
                .attr("transform", d => "translate(" + d.y + "," + d.x + ")");
            nodeUpdate.select('circle.node')
                .attr('r', 10)
                .style("fill", d => d._children ? "lightsteelblue" : "#fff")
                .attr('cursor', 'pointer');

            const nodeExit = node.exit().transition().duration(duration)
                .attr("transform", d => "translate(" + source.y + "," + source.x + ")")
                .remove();
            nodeExit.select('circle').attr('r', 1e-6);
            nodeExit.select('text').style('fill-opacity', 1e-6);

            const link = svg.selectAll('path.link').data(links, d => d.id);
            const linkEnter = link.enter().insert('path', "g")
                .attr("class", "link")
                .attr('d', d => {
                    const o = {x: source.x0, y: source.y0};
                    return diagonal(o, o);
                });

            const linkUpdate = linkEnter.merge(link);
            linkUpdate.transition().duration(duration).attr('d', d => diagonal(d, d.parent));
            link.exit().transition().duration(duration)
                .attr('d', d => {
                    const o = {x: source.x, y: source.y};
                    return diagonal(o, o);
                })
                .remove();

            let minX = 0, maxX = 0, minY = 0, maxY = 0;
            nodes.forEach(d => {
                minX = Math.min(minX, d.x);
                maxX = Math.max(maxX, d.x);
                minY = Math.min(minY, d.y);
                maxY = Math.max(maxY, d.y);
            });

            const newWidth = maxY - minY + margin.left + margin.right;
            const newHeight = maxX - minX + margin.top + margin.bottom;
            d3.select("svg").transition().duration(duration)
                .attr("width", newWidth)
                .attr("height", newHeight);
            svg.transition().duration(duration)
                .attr("transform", "translate(" + margin.left + "," + (-minX + margin.top) + ")");

            nodes.forEach(d => {
                d.x0 = d.x;
                d.y0 = d.y;
            });

            function diagonal(s, d) {
                return `M ${s.y} ${s.x}
                        C ${(s.y + d.y) / 2} ${s.x},
                          ${(s.y + d.y) / 2} ${d.x},
                          ${d.y} ${d.x}`
            }

            function click(d) {
                const detailsDiv = d3.select("#details");
                detailsDiv.html(`
                    <h3>${d.data.name}</h3>
                    <p><strong>Original Name:</strong> ${d.data.original_name}</p>
                    <p><strong>Leaf Snarls in Subtree:</strong> ${d.data.num_leaf_snarls}</p>
                `);
                if (d.children) {
                    d._children = d.children;
                    d.children = null;
                } else {
                    d.children = d._children;
                    d._children = null;
                }
                update(d);
            }
        }
    </script>
</body>
</html>)HTML";
    
    out.close();
}

int main(int argc, char* argv[]) {
    std::string json_path = "/private/groups/migalab/shnegi/vg_anchors_project/vg_anchors/snarl_tree_builder_cpp/snarl_tree_map.json";
    std::string output_path = "/private/groups/migalab/shnegi/vg_anchors_project/vg_anchors/snarl_tree_builder_cpp/snarl_tree.html";
    int initial_depth = 3;
    int filter_depth = 3;
    int min_leaf_snarls = 1;
    
    // Parse arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
            json_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            output_path = argv[++i];
        } else if ((arg == "-h" || arg == "--help")) {
            std::cout << "Usage: " << argv[0] << " [-i INPUT] [-o OUTPUT]\n";
            std::cout << "  -i, --input   Path to snarl tree JSON (default: snarl_tree_map.json)\n";
            std::cout << "  -o, --output  Path to output HTML (default: snarl_tree.html)\n";
            return 0;
        }
    }
    
    try {
        auto total_start = std::chrono::high_resolution_clock::now();
        
        // Load JSON
        SnarlTreeMap tree_map;
        tree_map.load_from_json(json_path);
        
        // Build tree
        TreeBuilder builder(tree_map);
        auto tree = builder.build_tree();
        
        // Serialize to JSON string
        std::cerr << "Serializing tree to JSON..." << std::endl;
        auto json_start = std::chrono::high_resolution_clock::now();
        
        std::ostringstream json_stream;
        builder.write_tree_json(json_stream, *tree);
        std::string tree_json = json_stream.str();
        
        auto json_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> json_time = json_end - json_start;
        std::cerr << "JSON serialization completed in " << std::fixed << std::setprecision(2)
                  << json_time.count() << " seconds" << std::endl;
        
        // Write HTML
        write_html_template(output_path, tree_json, initial_depth, filter_depth, min_leaf_snarls);
        
        auto total_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> total_time = total_end - total_start;
        
        std::cerr << "✓ Generated visualization at: " << output_path << std::endl;
        std::cerr << "Total time: " << std::fixed << std::setprecision(2) 
                  << total_time.count() << " seconds" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

