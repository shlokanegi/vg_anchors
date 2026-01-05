#!/usr/bin/env python3
"""Locate unique minimizers from walk_uniques.tsv inside a GFA assembly."""

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

Complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def parse_walk_tsv(path: Path) -> List[str]:
    """Return the ordered list of minimizer sequences."""
    minimizers: List[str] = []
    with path.open() as tsv:
        header = tsv.readline()
        if not header:
            return minimizers
        for idx, line in enumerate(tsv, start=2):
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4:
                raise ValueError(f"{path}:{idx}: expected 4 columns, got {len(parts)}")
            minimizer_seq = parts[0]
            minimizers.append(minimizer_seq)
    return minimizers


def parse_gfa_segments(path: Path) -> Dict[str, str]:
    segments: Dict[str, str] = {}
    with path.open() as gfa:
        for idx, line in enumerate(gfa, start=1):
            if not line or line[0] != "S":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                raise ValueError(f"{path}:{idx}: malformed S-line")
            segments[parts[1]] = parts[2]
    return segments


def parse_gfa_links(path: Path) -> List[Tuple[str, str, str, str]]:
    links: List[Tuple[str, str, str, str]] = []
    with path.open() as gfa:
        for idx, line in enumerate(gfa, start=1):
            if not line or line[0] != "L":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                raise ValueError(f"{path}:{idx}: malformed L-line")
            links.append((parts[1], parts[2], parts[3], parts[4]))
    return links


def revcomp(seq: str) -> str:
    return seq.translate(Complement)[::-1]


def build_minimizer_index(minimizers: List[str]) -> Dict[int, Dict[str, int]]:
    by_length: Dict[int, Dict[str, int]] = defaultdict(dict)
    for idx, seq in enumerate(minimizers):
        by_length[len(seq)][seq] = idx
    return by_length


def record_hit(summary, contig: str, idx: int, offset: int, length: int):
    if contig not in summary:
        summary[contig] = {"indices": [], "offsets": [], "lengths": [], "seen": set()}
    entry = summary[contig]
    if idx in entry["seen"]:
        return
    entry["seen"].add(idx)
    entry["indices"].append(str(idx))
    entry["offsets"].append(str(offset))
    entry["lengths"].append(str(length))


def scan_within_contigs(
    segments: Dict[str, str],
    min_index: Dict[int, Dict[str, int]],
    summary,
):
    for contig, seq in segments.items():
        seq_len = len(seq)
        for length, lookup in min_index.items():
            if seq_len < length:
                continue
            for pos in range(seq_len - length + 1):
                substring = seq[pos : pos + length]
                idx = lookup.get(substring)
                if idx is None:
                    continue
                record_hit(summary, contig, idx, pos, length)


def oriented_offset(offset: int, match_len: int, orient: str, seq_len: int) -> int:
    if orient == "+":
        return offset
    return seq_len - offset - match_len


def scan_links(
    segments: Dict[str, str],
    links: List[Tuple[str, str, str, str]],
    min_index: Dict[int, Dict[str, int]],
    summary,
):
    for src, src_orient, dst, dst_orient in links:
        src_seq = segments.get(src)
        dst_seq = segments.get(dst)
        if src_seq is None or dst_seq is None:
            continue
        src_seq_oriented = src_seq if src_orient == "+" else revcomp(src_seq)
        dst_seq_oriented = dst_seq if dst_orient == "+" else revcomp(dst_seq)
        len_src = len(src_seq)
        len_dst = len(dst_seq)

        for length, lookup in min_index.items():
            if length <= 1:
                continue
            max_split = min(length - 1, len(src_seq_oriented))
            for split in range(1, max_split + 1):
                remain = length - split
                if len(dst_seq_oriented) < remain:
                    continue
                candidate = (
                    src_seq_oriented[-split:] + dst_seq_oriented[:remain]
                )
                idx = lookup.get(candidate)
                if idx is None:
                    continue
                src_off_oriented = len(src_seq_oriented) - split
                dst_off_oriented = 0
                src_offset = oriented_offset(src_off_oriented, split, src_orient, len_src)
                dst_offset = oriented_offset(dst_off_oriented, remain, dst_orient, len_dst)
                record_hit(summary, src, idx, src_offset, split)
                record_hit(summary, dst, idx, dst_offset, remain)


def summarize(
    minimizers: List[str],
    segments: Dict[str, str],
    links: List[Tuple[str, str, str, str]],
):
    summary = {}
    min_index = build_minimizer_index(minimizers)
    scan_within_contigs(segments, min_index, summary)
    scan_links(segments, links, min_index, summary)
    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Map graph-unique minimizers to contig sequences"
    )
    parser.add_argument("--walk-tsv", required=True, type=Path, help="walk_uniques.tsv")
    parser.add_argument("--gfa", required=True, type=Path, help="assembly GFA")
    parser.add_argument("--output", required=True, type=Path, help="output TSV path")
    args = parser.parse_args()

    minimizers = parse_walk_tsv(args.walk_tsv)
    if not minimizers:
        raise SystemExit("No minimizers found in walk TSV")
    segments = parse_gfa_segments(args.gfa)
    if not segments:
        raise SystemExit("No segments present in GFA")
    links = parse_gfa_links(args.gfa)

    summary = summarize(minimizers, segments, links)
    with args.output.open("w") as out:
        out.write(
            "contig\tnumber_of_graph_unique_minimizers\tminimizer_indices\tcontig_offsets\tminimizer_lengths\n"
        )
        for contig, record in summary.items():
            out.write(
                f"{contig}\t{len(record['seen'])}\t{','.join(record['indices'])}\t"
                f"{','.join(record['offsets'])}\t{','.join(record['lengths'])}\n"
            )


if __name__ == "__main__":
    main()

