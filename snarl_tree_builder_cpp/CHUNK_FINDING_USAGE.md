# Chunk Point Finding Feature

## Overview

The `chunk_point_finding` tool supports automatic generation of chunk points for parallel subgraph extraction. It intelligently subdivides the variation graph into chunks containing approximately a target number of leaf snarls (default: 100,000).

**New in this version**: Fast mode with JSON loading to iterate quickly on different chunk sizes!

## Three Operating Modes

### Mode 1: Build Snarl Tree Only
Build and save the tree for later use:
```bash
./chunk_point_finding -g graph.pg -i index.dist -o tree.json -n 64
```

### Mode 2: Fast Chunk Generation ⭐ Recommended, but you should have the snarl tree already built
Load existing tree and generate chunks (**saves ~100 seconds**):
```bash
./chunk_point_finding -j tree.json -g graph.pg -i index.dist -c chunks.tsv -t 100000
```

### Mode 3: Build Tree + Generate Chunks
Do everything in one step:
```bash
./chunk_point_finding -g graph.pg -i index.dist -c chunks.tsv -n 64
```

## Algorithm

The chunking algorithm works as follows:

1. **Start with root-level chains**: Process all top-level chains (typically ~80)

2. **Skip empty chains**: Chains with 0 leaf snarls are excluded

3. **For each non-empty chain**:
   - If total leaf snarls ≤ target: Add the entire chain as a single chunk
   - If total leaf snarls > target: Subdivide the chain

4. **Subdivision process**:
   - Iterate through children in order (preserving snarl decomposition order)
   - Accumulate leaf snarl counts from both snarls and chains
   - When cumulative count ≥ target OR at last child:
     - Create chunk from last boundary to current position
     - Reset counter and continue (if not last child)

5. **Output**: TSV file with columns (start_node, end_node, leaf_snarls)

## Command-Line Options

### Required (varies by mode)
- **Mode 1**: `-g`, `-i`
- **Mode 2**: `-j`, `-g`, `-i`, `-c`
- **Mode 3**: `-g`, `-i`, `-c`

### All Options
- `-g, --graph`: Path to variation graph (.pg file)
- `-i, --index`: Path to snarl distance index (.dist file)
- `-j, --input-json`: Load existing snarl tree JSON (enables fast mode)
- `-o, --output-json`: Path to output snarl tree JSON (default: snarl_tree_map.json)
- `-n, --num-threads`: Number of threads to use (default: auto-detect)
- `-c, --chunk-points`: Path to output chunk points TSV file
- `-t, --target-chunk-size`: Target leaf snarls per chunk (default: 100000)

## Output Format

The output is a TSV (tab-separated values) file with three columns:

```
start_node	end_node	leaf_snarls
12345	23456	98543
23456	34567	105234
34567	45678	87456
```

Each row represents a chunk:
- **start_node**: Starting boundary node ID
- **end_node**: Ending boundary node ID
- **leaf_snarls**: Number of leaf snarls in this chunk (useful for validation)

## Example Output

### First Run (Full Build)
```
Loading graph: 45.65 seconds
Loading index: 8.70 seconds
Running with 64 parallel threads.
Building snarl tree: 136.54 seconds (produced 9720133 entries)
Processing 80 root-level chains for chunking...
Generated 116 chunk points (target: ~100000 leaf snarls per chunk)
Chunk points written to: chunk_points.tsv
```

### Fast Mode (Load Existing Tree)
```
Loading graph: 47.18 seconds
Loading index: 8.92 seconds
Loaded snarl tree from JSON: 26.49 seconds (9720133 entries)
Processing 80 root-level chains for chunking...
Generated 172 chunk points (target: ~50000 leaf snarls per chunk)
Chunk points written to: chunks_50k.tsv
```

**Time saved**: ~100 seconds by loading JSON instead of rebuilding!

## Use Case: Parallel Subgraph Extraction

### Bash: Simple iteration
```bash
while IFS=$'\t' read -r start end leaves; do
    if [ "$start" != "start_node" ]; then  # Skip header
        echo "Processing chunk: $start to $end ($leaves leaf snarls)"
        # vg find -x graph.xg -n $start:$end > chunk_${start}_${end}.vg
    fi
done < chunk_points.tsv
```

### Python: With validation
```python
import pandas as pd
from multiprocessing import Pool

chunks = pd.read_csv('chunk_points.tsv', sep='\t')

# Validate chunks
print(f"Total chunks: {len(chunks)}")
print(f"Total leaf snarls: {chunks['leaf_snarls'].sum()}")
print(f"Average per chunk: {chunks['leaf_snarls'].mean():.0f}")
print(f"Min/Max: {chunks['leaf_snarls'].min()} / {chunks['leaf_snarls'].max()}")

def process_chunk(row):
    start, end, leaf_count = row['start_node'], row['end_node'], row['leaf_snarls']
    print(f"Processing {leaf_count} leaf snarls in chunk {start}→{end}")
    # Your subgraph extraction code here
    
with Pool(16) as pool:
    pool.map(process_chunk, [row for _, row in chunks.iterrows()])
```

## Workflow: Iterating on Chunk Sizes

1. **Build tree once** (saves time later):
   ```bash
   ./chunk_point_finding -g graph.pg -i index.dist -o tree.json -n 64
   ```

2. **Try different chunk sizes** (fast iterations):
   ```bash
   # Large chunks (fewer, more work each)
   ./chunk_point_finding -j tree.json -g graph.pg -i index.dist -c chunks_200k.tsv -t 200000
   
   # Small chunks (many, less work each)
   ./chunk_point_finding -j tree.json -g graph.pg -i index.dist -c chunks_50k.tsv -t 50000
   
   # Default
   ./chunk_point_finding -j tree.json -g graph.pg -i index.dist -c chunks_100k.tsv -t 100000
   ```

3. **Validate results**:
   ```bash
   for file in chunks_*.tsv; do
       echo "$file:"
       awk 'NR>1 {sum+=$3; count++} END {print "  Chunks:", count, "| Avg leaves:", int(sum/count)}' $file
   done
   ```

## Benefits

1. **Fast iteration**: Load JSON once, experiment with chunk sizes quickly
2. **Automatic**: No manual selection of chunk boundaries
3. **Balanced**: Each chunk contains approximately the target leaf snarls
4. **Ordered**: Chunks follow the snarl decomposition order
5. **Validated**: Includes leaf snarl count for each chunk
6. **Efficient**: Generated in-memory, minimal overhead

## Performance Tips

- **Use Mode 2** when iterating on chunk sizes (saves ~100 seconds per run)
- **Larger chunks** (e.g., 200k): Fewer chunks, less overhead, more memory per chunk
- **Smaller chunks** (e.g., 50k): More chunks, better load balancing, less memory per chunk
- **Default (100k)**: Good balance for most use cases
- **Validate first**: Check the `leaf_snarls` column before running heavy downstream workflows

## Implementation Details

- Chunks are created by traversing children in order
- Both snarls and chains can be chunk boundaries
- Empty chains (0 leaf snarls) are automatically skipped
- Last chunk in subdivided chains may be smaller than target
- Chains already < target size are kept as single chunks
