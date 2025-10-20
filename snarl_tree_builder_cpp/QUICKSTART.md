# Quick Start Guide

## Choose Your Tool

- **`construct_snarl_tree_map`**: Basic snarl tree builder
- **`chunk_point_finding`**: Advanced tool with chunking and fast mode ⭐ **Recommended**

## Build in 3 Steps

1. **Navigate to the directory:**
   ```bash
   cd vg_anchors/snarl_tree_builder_cpp
   ```

2. **Build the executables:**
   ```bash
   make
   ```

3. **Run chunk_point_finding (recommended):**
   ```bash
   ./chunk_point_finding -g <graph.pg.vg> -i <graph.pg.dist> -n 64
   ```

## Example: Basic Snarl Tree

```bash
# Build
make

# Run on your data
./construct_snarl_tree_map \
    -g /path/to/your/graph.pg.vg \
    -i /path/to/your/index.pg.dist \
    -o snarl_tree.json

# Output:
# Loading graph: 0.39 seconds
# Loading index: 0.04 seconds
# Running with 32 parallel threads.
# Building snarl tree: 0.84 seconds (produced 66704 entries)
```

## Example: Chunk Points (Recommended Workflow)

### First Run: Build Tree + Generate Chunks
```bash
./chunk_point_finding \
    -g /path/to/graph.pg.vg \
    -i /path/to/index.pg.dist \
    -c chunk_points.tsv \
    -n 64 \
    -t 100000

# Output:
# Loading graph: 45.65 seconds
# Loading index: 8.70 seconds
# Running with 64 parallel threads.
# Building snarl tree: 136.54 seconds (produced 9720133 entries)
# Processing 80 root-level chains for chunking...
# Generated 116 chunk points (target: ~100000 leaf snarls per chunk)
# Chunk points written to: chunk_points.tsv
```

### Subsequent Runs: Fast Mode (Load Existing Tree)
```bash
./chunk_point_finding \
    -j snarl_tree_map.json \
    -g /path/to/graph.pg.vg \
    -i /path/to/index.pg.dist \
    -c chunk_points_50k.tsv \
    -t 50000

# Output:
# Loading graph: 47.18 seconds
# Loading index: 8.92 seconds
# Loaded snarl tree from JSON: 26.49 seconds (9720133 entries)
# Processing 80 root-level chains for chunking...
# Generated 172 chunk points (target: ~50000 leaf snarls per chunk)
# Chunk points written to: chunk_points_50k.tsv
#
# Time saved: ~100 seconds by loading JSON instead of rebuilding!
```

## Expected Output Formats

### Snarl Tree JSON
```json
{
    "s80156932-80156928": [
        ["c80156930-80156929", "80156928"],
        1
    ],
    ...
}
```

### Chunk Points TSV
```
start_node	end_node	leaf_snarls
115630569	116793163	100000
116793163	125555905	100000
125555905	119661771	100000
```

## Troubleshooting

**Build fails?**
- Check that `../libbdsg/` exists and is built
- Try: `cd ../libbdsg && make`

**Segmentation fault?**
- Verify your graph and index files are not corrupted
- Ensure graph and index files match (same graph)

## Visualizing the Tree

After generating the JSON, create an interactive HTML visualization:

```bash
# Generates snarl_tree.html with D3.js visualization
./visualize_snarl_tree
```

Opens the HTML file in a browser to:
- Explore the tree interactively
- Filter by depth and leaf snarl count
- Click nodes for details

For more details, see [README.md](README.md)

