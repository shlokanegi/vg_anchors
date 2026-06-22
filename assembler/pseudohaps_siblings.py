#!/usr/bin/env python3
"""
pseudohaps_siblings.py

Part 1 of the anchor-based pseudohaplotype pipeline: compute the anchor/snarl
sibling relationships between assembled contigs and write them to
``{output_prefix}_sibling_contigs.tsv``.

This is the expensive half — it loads each chunk's Shasta AssemblyDetails /
AnchorsFromJson and the extended snarl TSV to build per-contig snarl profiles.
Run it once; ``pseudohaps_build.py`` then consumes the sibling TSV to construct
hap1/hap2 (and can be re-run cheaply while iterating on the walk logic).

See pseudohaps_anchors.py for the shared implementation.
"""

from __future__ import annotations

import argparse
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set

import pseudohaps as ph
import pseudohaps_anchors as pa


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gfa", help="Assembled GFA (chunk-prefixed contig names)")
    parser.add_argument("-c", "--chunk-points", required=True)
    parser.add_argument("-o", "--output-prefix", required=True)
    parser.add_argument(
        "--region-dir",
        default=None,
        help="Region root with per-chunk subdirs (0/, 1/, ...). "
             "Default: inferred as parent of combined/ from the GFA path.",
    )
    parser.add_argument(
        "--chunk",
        action="append",
        default=[],
        metavar="ID=SHASTA_RUN_DIR",
        help="Optional manual Shasta run dir override (repeat per chunk).",
    )
    parser.add_argument(
        "--extended-tsv",
        action="append",
        default=[],
        metavar="ID=PATH",
        help="Optional override for subgraph.anchors.json.subgraph.sizes.extended.tsv",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Snakemake YAML config; SKIP_CHROMOSOMES are excluded from anchor loading.",
    )
    parser.add_argument(
        "--skip-chromosome",
        action="append",
        default=[],
        help="Additional chromosome to skip (repeatable). Merged with SKIP_CHROMOSOMES from --config.",
    )
    parser.add_argument("--min-match-fraction", type=float, default=0.5)
    parser.add_argument("--min-shared-snarls", type=int, default=2)
    parser.add_argument("--min-snarl-fraction", type=float, default=0.5)
    args = parser.parse_args()

    gfa_path = Path(args.gfa)
    region_dir = pa.infer_region_dir(
        gfa_path,
        Path(args.region_dir) if args.region_dir else None,
    )
    shasta_suffix = pa.extract_shasta_run_suffix(gfa_path)
    chunk_to_chrom_all = ph.parse_chunk_points_file(args.chunk_points)

    skip_chromosomes: Set[str] = set(args.skip_chromosome)
    if args.config:
        skip_chromosomes.update(pa.load_skip_chromosomes_from_config(Path(args.config)))

    chunk_to_chrom, skipped_chunks = pa.apply_skip_chromosomes(
        chunk_to_chrom_all, skip_chromosomes
    )
    chunk_ids = sorted(chunk_to_chrom.keys(), key=lambda x: int(x))

    extended_overrides: Dict[str, Path] = {}
    for spec in args.extended_tsv:
        cid, p = spec.split("=", 1)
        extended_overrides[cid.strip()] = Path(p.strip())

    print(f"\nRegion directory: {region_dir}")
    if shasta_suffix:
        print(f"Shasta run suffix from GFA name: {shasta_suffix}")
    if skip_chromosomes:
        print(f"Skipping chromosomes: {sorted(skip_chromosomes)}")
    if skipped_chunks:
        by_chrom: Dict[str, list] = defaultdict(list)
        for cid, chrom in skipped_chunks:
            by_chrom[chrom].append(cid)
        for chrom in sorted(by_chrom):
            ids = sorted(by_chrom[chrom], key=lambda x: int(x))
            print(f"  Excluded {len(ids)} chunk(s) on {chrom}: {ids[0]}..{ids[-1]}")
    if not chunk_ids:
        raise ValueError(
            "No chunks remain after applying SKIP_CHROMOSOMES filters. "
            "Check --config and --skip-chromosome settings."
        )
    print(f"Chunk IDs for anchor loading: {len(chunk_ids)} ({chunk_ids[0]}..{chunk_ids[-1]})")

    print("\nLoading chunk anchor data...")
    t0 = time.time()
    chunk_data: Dict[str, pa.ChunkAnchorData] = {}

    if args.chunk:
        chunk_specs = pa.parse_chunk_specs(args.chunk)
        for cid, shasta_dir in chunk_specs:
            ext = extended_overrides.get(cid) or pa.discover_extended_tsv(region_dir, cid)
            chunk_data[cid] = pa.load_chunk_anchor_data(cid, shasta_dir, ext)
            print(f"  chunk {cid}: {len(chunk_data[cid].profiles)} contig profiles from {shasta_dir}")
    else:
        discovered = pa.discover_chunks_from_region(region_dir, chunk_ids, shasta_suffix)
        for cid, shasta_dir, ext_tsv in discovered:
            ext = extended_overrides.get(cid, ext_tsv)
            chunk_data[cid] = pa.load_chunk_anchor_data(cid, shasta_dir, ext)
            print(f"  chunk {cid}: {len(chunk_data[cid].profiles)} contig profiles from {shasta_dir}")
    print(f"  Done in {time.time() - t0:.2f}s")

    print("\nLoading GFA node set...")
    t0 = time.time()
    graph = ph.load_gfa(args.gfa)
    print(f"  {len(graph.all_nodes)} nodes in {time.time() - t0:.2f}s")

    # Profile fusion contigs (stitched cross-chunk joins, names with '+') from
    # their component contigs so they participate in sibling detection too.
    fusion_added = pa.attach_fusion_profiles(chunk_data, graph.all_nodes)
    print(f"  profiled {len(fusion_added)} fusion contigs from their components")

    # All GFA nodes that have anchor profiles (siblings live within the assembly).
    anchor_gfa_nodes: Set[str] = set()
    for cdata in chunk_data.values():
        anchor_gfa_nodes.update(cdata.profiles.keys())
    anchor_gfa_nodes &= graph.all_nodes

    print("\nBuilding anchor-based sibling index...")
    t0 = time.time()
    sibling_dict, pair_records = pa.build_sibling_index(
        chunk_data,
        anchor_gfa_nodes,
        graph,
        min_match_fraction=args.min_match_fraction,
        min_shared_snarls=args.min_shared_snarls,
        min_snarl_fraction=args.min_snarl_fraction,
    )
    sibling_tsv = f"{args.output_prefix}_sibling_contigs.tsv"
    pa.write_sibling_pairs_tsv(sibling_tsv, pair_records)
    print(f"  {len(pair_records)} undirected sibling pairs")
    print(f"  {sum(1 for v in sibling_dict.values() if v)} contigs with >=1 sibling")
    print(f"  Wrote {sibling_tsv} in {time.time() - t0:.2f}s")

    print("\nDone.")


if __name__ == "__main__":
    main()
