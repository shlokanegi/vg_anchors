# Snarl Tree Builder (C++)

A high-performance C++ implementation for constructing hierarchical snarl tree maps from variation graphs. This tool uses multithreading to efficiently process large graphs by traversing snarl decompositions in parallel.

## Overview

This program:
- Loads a variation graph and its snarl distance index
- Traverses the snarl decomposition to build a hierarchical tree
- Counts leaf snarls in each subtree
- Uses C++ multithreading for parallel processing at depth 2 of the tree
- Outputs a JSON file representing the complete snarl tree structure

## Features

- **Multithreaded**: Uses C++ threads (not multiprocessing) for efficient parallel execution
- **No Serialization Issues**: Passes handles directly between threads without pickling
- **High Performance**: C++ implementation with O3 optimization
- **Memory Efficient**: Threads share memory space, reducing overhead
- **Thread Pool Pattern**: Distributes work efficiently across available CPU cores

## Requirements

- **C++ Compiler**: g++ with C++17 support (or equivalent)
- **libbdsg**: The library must be installed in `../libbdsg/`
- **CMake**: Version 3.10+ (optional, if using CMake build)

## Building

### Option 1: Using Make (Recommended)

```bash
make
```

This will compile the executable `construct_snarl_tree_map` directly.

### Option 2: Using CMake

```bash
mkdir build
cd build
cmake ..
make
```

The executable will be created in the `build/` directory.

### Clean Build

```bash
make clean
```

## Usage

```bash
./construct_snarl_tree_map -g <graph_file> -i <index_file> [-o <output_json>]
```

### Arguments

- `-g, --graph`: Path to the variation graph file (.pg or .vg file)
- `-i, --index`: Path to the snarl distance index (.dist file)
- `-o, --output-json`: Path to the output JSON file (default: `snarl_tree_map.json`)
- `-h, --help`: Show help message

### Example

```bash
./construct_snarl_tree_map \
    -g /path/to/graph.pg.vg \
    -i /path/to/graph.pg.dist \
    -o snarl_tree_output.json
```

## Output Format

The output is a JSON file with the following structure:

```json
{
    "snarl_or_chain_id": [
        ["child_id_1", "child_id_2", ...],
        total_leaf_snarls_count
    ],
    ...
}
```

Where:
- **Keys**: Snarl/chain IDs (format: `s{node1}-{node2}` for snarls, `c{node1}-{node2}` for chains)
- **Values**: Array containing:
  - List of children IDs
  - Total count of leaf snarls in the subtree

### Example Entry

```json
"s80156932-80156928": [
    ["c80156930-80156929", "80156928"],
    1
]
```

This represents a snarl with boundaries at nodes 80156932 and 80156928, having 2 children, and containing 1 leaf snarl in its subtree.

## Performance

### Parallelization Strategy

The tool uses a depth-based parallelization approach:

1. Traverses the tree to depth 2
2. Collects all snarls/chains at that depth
3. Creates worker threads (up to CPU count)
4. Each thread processes a subset of subtrees
5. Results are merged thread-safely
6. Top levels of the tree are reconstructed

### Threading vs Multiprocessing

This C++ implementation uses **threading** instead of Python's multiprocessing:
- ✅ **No serialization overhead**: Handles are passed directly
- ✅ **Shared memory**: Graph and index loaded once
- ✅ **Lower overhead**: Thread creation is faster than process forking
- ✅ **Better cache locality**: Threads share CPU cache

### Typical Performance

For a 30Mb chromosome graph (~65k handles at depth 2):
- **Compilation**: ~5 seconds
- **Total runtime**: ~3 seconds (parallelized across all available cores)

### Output

The tool provides minimal, clean logging:
```
Loading graph: 0.39 seconds
Loading index: 0.04 seconds
Running with 32 parallel threads.
Building snarl tree: 0.84 seconds (produced 66704 entries)
```

## Implementation Details

### Key Classes and Functions

- **`SnarlTreeBuilder`**: Main class managing the decomposition
  - `load()`: Deserializes graph and index
  - `build_parallel()`: Orchestrates parallel tree construction
  - `process_subtree_worker()`: Worker function for each thread
  - `collect_handles_at_depth()`: Gathers handles for parallelization
  - `build_top_tree()`: Reconstructs upper tree levels
  - `traverse_decomposition()`: Recursive tree traversal
  - `is_leaf_snarl()`: Determines if a snarl is a leaf

### Thread Safety

- Each worker thread builds its own local map
- A mutex protects the merge operation into the shared map
- No race conditions or data corruption

### Comparison to Python Version

| Aspect | Python | C++ |
|--------|--------|-----|
| Speed | Slower | ~10x faster |
| Memory | Higher (process copies) | Lower (shared memory) |
| Handle passing | Pickle serialization required | Direct passing |
| Dependencies | Python, bdsg-python | C++17, libbdsg |
| Ease of use | Easier syntax | More verbose |

## Troubleshooting

### Compilation Errors

**Error: Cannot find bdsg headers**
```
Solution: Ensure libbdsg is installed in ../libbdsg/
Check paths in Makefile or CMakeLists.txt
```

**Error: Undefined references during linking**
```
Solution: Verify all libbdsg libraries are built
Check that library paths are correct in LDFLAGS
```

### Runtime Errors

**Error: Cannot open graph/index file**
```
Solution: Check file paths are correct and files exist
Ensure you have read permissions
```

**Error: Segmentation fault**
```
Solution: This may indicate corrupted input files
Try with a smaller test graph first
Check that graph and index match
```

## Development

### Code Structure

```
snarl_tree_builder_cpp/
├── construct_snarl_tree_map.cpp  # Main implementation
├── CMakeLists.txt                 # CMake build configuration
├── Makefile                       # Direct make configuration
└── README.md                      # This file
```

### Modifying Parallelization Depth

To change the depth at which parallelization occurs, modify the constant:

```cpp
const int PARALLELIZATION_DEPTH = 2;  // Line 32 in .cpp file
```

Higher values = more granular parallelization but more overhead.

## License

This tool is part of the vg_anchors project.

## Chunk Point Finding

The `chunk_point_finding` tool automatically generates chunk points for parallel subgraph extraction. It supports three operating modes for maximum flexibility.

### Quick Start

**Fast Mode (recommended for iterations):**
```bash
# First time: build the snarl tree (saves to snarl_tree_map.json)
./chunk_point_finding -g graph.pg -i index.dist -n 64

# Subsequent runs: load existing tree (much faster!)
./chunk_point_finding -j snarl_tree_map.json -g graph.pg -i index.dist -c chunks.tsv -t 100000
```

### Three Operating Modes

#### Mode 1: Build Snarl Tree Only
Build and save the snarl tree for later use:
```bash
./chunk_point_finding -g graph.pg -i index.dist -o tree.json -n 64
```

#### Mode 2: Fast Chunk Generation (Recommended)
Load existing snarl tree and generate chunks (**saves ~100 seconds**):
```bash
./chunk_point_finding -j tree.json -g graph.pg -i index.dist -c chunks.tsv
```

#### Mode 3: Build Tree + Generate Chunks
Do everything in one step:
```bash
./chunk_point_finding -g graph.pg -i index.dist -c chunks.tsv -n 64
```

### Features

- **Three flexible modes**: Build tree, load tree, or do both
- **Fast iteration**: Reuse existing tree when experimenting with chunk sizes
- **Automatic chunking**: Intelligently subdivides the graph into balanced chunks
- **Configurable size**: Specify target number of leaf snarls per chunk (default: 100k)
- **Order-preserving**: Chunks follow the snarl decomposition order
- **Skip empty chains**: Automatically excludes chains with 0 leaf snarls

### Command-Line Arguments

- `-g, --graph`: Path to variation graph (.pg file) - **required**
- `-i, --index`: Path to snarl distance index (.dist file) - **required**
- `-j, --input-json`: Load existing snarl tree JSON (enables fast mode)
- `-o, --output-json`: Path to output snarl tree JSON (default: snarl_tree_map.json)
- `-n, --num-threads`: Number of threads to use (default: auto-detect)
- `-c, --chunk-points`: Path to output chunk points TSV file
- `-t, --target-chunk-size`: Target leaf snarls per chunk (default: 100000)

### Output Format

TSV file with three columns:
```
start_node	end_node	leaf_snarls
12345	23456	98543
23456	34567	105234
34567	45678	87456
```

- **start_node**: Starting boundary node ID of the chunk
- **end_node**: Ending boundary node ID of the chunk
- **leaf_snarls**: Number of leaf snarls in this chunk

### Performance

For a large human chromosome graph (~9.7M entries):
- **Full build**: ~3 minutes (54s load + 130s build + chunking)
- **Fast mode**: ~2.5 minutes (56s load + 26s JSON load + chunking)
- **Time saved**: ~100 seconds per iteration

### Algorithm

1. Iterate through all root-level chains (typically ~80)
2. Skip chains with 0 leaf snarls (no useful variation)
3. For chains with ≤ target leaf snarls: keep as single chunk
4. For larger chains: subdivide by traversing children in order
   - Accumulate leaf snarls from snarls and chains
   - Create chunk boundary when target is reached
   - Handle both snarls and chains as potential boundaries
5. Output chunks with accurate leaf snarl counts

### Example Use Cases

**Python: Parallel subgraph extraction**
```python
import pandas as pd
from multiprocessing import Pool

chunks = pd.read_csv('chunk_points.tsv', sep='\t')

def process_chunk(row):
    start, end, leaf_count = row['start_node'], row['end_node'], row['leaf_snarls']
    print(f"Processing chunk with {leaf_count} leaf snarls")
    # Your subgraph extraction code here
    
with Pool(16) as pool:
    pool.map(process_chunk, [row for _, row in chunks.iterrows()])
```

**Bash: Simple iteration**
```bash
while IFS=$'\t' read -r start end leaves; do
    if [ "$start" != "start_node" ]; then  # Skip header
        echo "Processing chunk: $start to $end ($leaves leaf snarls)"
        # vg find -x graph.xg -n $start:$end > chunk_${start}_${end}.vg
    fi
done < chunk_points.tsv
```

### Tips for Iterating

1. **Build the tree once**: Use Mode 1 to create `tree.json`
2. **Experiment with chunk sizes**: Use Mode 2 with different `-t` values
3. **Validate chunks**: Check the `leaf_snarls` column to ensure balanced distribution
4. **Use threads**: Specify `-n 64` or higher for large graphs

## Visualization Tool

The package also includes a fast C++ visualization tool:

```bash
./visualize_snarl_tree -i snarl_tree_map.json -o snarl_tree.html
```

This generates an interactive HTML visualization with D3.js. Performance for large trees:
- Loads 9.7M entries in ~16 seconds
- Builds tree structure in ~11 seconds  
- Total time: ~56 seconds

The visualization includes:
- Interactive expand/collapse nodes
- Filter by depth and minimum leaf snarls
- Click nodes to view details
- Automatic layout adjustment

## Contact

For issues or questions, please refer to the main vg_anchors project documentation.

