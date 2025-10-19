# Quick Start Guide

## Build in 3 Steps

1. **Navigate to the directory:**
   ```bash
   cd vg_anchors/snarl_tree_builder_cpp
   ```

2. **Build the executable:**
   ```bash
   make
   ```

3. **Run it:**
   ```bash
   ./construct_snarl_tree_map -g <graph.pg.vg> -i <graph.pg.dist> -o output.json
   ```

## Example

```bash
# Build
make

# Run on your data
./construct_snarl_tree_map \
    -g /path/to/your/graph.pg.vg \
    -i /path/to/your/index.pg.dist \
    -o snarl_tree.json

# The output will show:
# - Time to load graph and index
# - Number of threads used
# - Time to build snarl tree and entry count
```

## Expected Output

```
Loading graph: 0.39 seconds
Loading index: 0.04 seconds
Running with 32 parallel threads.
Building snarl tree: 0.84 seconds (produced 66704 entries)
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

