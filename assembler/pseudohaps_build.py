#!/usr/bin/env python3
"""
pseudohaps_build.py

Part 2 of the anchor-based pseudohaplotype pipeline: build hap1/hap2 from the
assembled GFA using a precomputed sibling TSV (from ``pseudohaps_siblings.py``).

This step does NOT touch the Shasta anchor data — it only needs the GFA, the
chunk-points TSV (to group paths by chromosome) and the sibling relationships,
so it is cheap to re-run while iterating on the hap2 walk.

  hap1: longest bp path through each assembled chain.
  hap2: sibling-aware graph walk (see pseudohaps_anchors.construct_hap2_path).

Outputs (prefix ``{output_prefix}``):
  _pseudohaplotypes.tsv, _hap1.fasta, _hap2.fasta, _ploidy_assignments.tsv
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Set

import pseudohaps as ph
import pseudohaps_anchors as pa


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gfa", help="Assembled GFA (chunk-prefixed contig names)")
    parser.add_argument("-c", "--chunk-points", required=True)
    parser.add_argument("-o", "--output-prefix", required=True)
    parser.add_argument(
        "-s", "--sibling-tsv", required=True,
        help="Sibling TSV produced by pseudohaps_siblings.py "
             "({prefix}_sibling_contigs.tsv).",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Snakemake YAML config; SKIP_CHROMOSOMES are excluded from chain grouping.",
    )
    parser.add_argument(
        "--skip-chromosome",
        action="append",
        default=[],
        help="Additional chromosome to skip (repeatable). Merged with SKIP_CHROMOSOMES from --config.",
    )
    args = parser.parse_args()

    chunk_to_chrom_all = ph.parse_chunk_points_file(args.chunk_points)

    skip_chromosomes: Set[str] = set(args.skip_chromosome)
    if args.config:
        skip_chromosomes.update(pa.load_skip_chromosomes_from_config(Path(args.config)))

    chunk_to_chrom, skipped_chunks = pa.apply_skip_chromosomes(
        chunk_to_chrom_all, skip_chromosomes
    )
    if skip_chromosomes:
        print(f"Skipping chromosomes: {sorted(skip_chromosomes)}")

    print("\nLoading GFA...")
    t0 = time.time()
    graph = ph.load_gfa(args.gfa)
    print(f"  {len(graph.all_nodes)} nodes in {time.time() - t0:.2f}s")

    print(f"\nLoading sibling relationships from {args.sibling_tsv}...")
    sibling_dict = pa.load_sibling_dict_from_tsv(args.sibling_tsv)
    print(f"  {sum(1 for v in sibling_dict.values() if v)} contigs with >=1 sibling")

    print("\nFinding paths and assembled chains...")
    t0 = time.time()
    paths = list(
        graph.find_all_start_to_end_paths(
            graph.get_potential_start_nodes(),
            graph.get_potential_end_nodes(),
        ).items()
    )
    chain_candidates = ph.group_paths_by_chromosome_and_assembled_chains(
        paths, chunk_to_chrom
    )
    assembled_chains = ph.select_longest_paths_per_assembled_chain_candidate(
        graph, chain_candidates
    )
    print(f"  {sum(len(v) for v in assembled_chains.values())} chains in {time.time() - t0:.2f}s")

    print("\nConstructing pseudohaplotypes (hap1=longest bp, hap2=sibling walk)...")
    t0 = time.time()
    pseudohaplotypes = pa.construct_pseudohaplotypes_anchors(
        graph, assembled_chains, sibling_dict
    )
    n_pairs = sum(len(v) for v in pseudohaplotypes.values())
    print(f"  {n_pairs} haplotype pairs in {time.time() - t0:.2f}s")

    print("\nWriting pseudohaplotypes TSV...")
    pa.write_pseudohaplotypes_tsv(
        f"{args.output_prefix}_pseudohaplotypes.tsv", pseudohaplotypes, graph
    )

    print("\nWriting FASTA files...")
    hap1_fasta, hap2_fasta = pa.write_haplotype_fastas(
        args.output_prefix, pseudohaplotypes, graph
    )
    print(f"  {hap1_fasta}, {hap2_fasta}")

    print("\nWriting ploidy assignments...")
    ploidy_rows = pa.compute_ploidy_assignments(pseudohaplotypes, graph, sibling_dict)
    with open(f"{args.output_prefix}_ploidy_assignments.tsv", "w") as out:
        out.write("chromosome\tassembled_chain_id\tcontig_id\tploidy\tcontext\tlength\n")
        for row in ploidy_rows:
            out.write("\t".join(str(x) for x in row) + "\n")

    print("\nWriting Bandage CSVs...")
    b1, b2 = pa.write_bandage_csvs(args.output_prefix, pseudohaplotypes, graph, sibling_dict)
    print(f"  {b1}, {b2}")

    print("\nDone.")


if __name__ == "__main__":
    main()
