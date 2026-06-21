#!/usr/bin/env python3
"""
pseudohaps_anchors.py

Build pseudohaplotypes from an assembled GFA using anchor/snarl-based sibling
contig detection (Shasta AssemblyDetails + AnchorsFromJson + extended snarl TSV).

Chunk Shasta directories are auto-discovered from:
  - chunk IDs in the chunk-points TSV (row index 0, 1, 2, ... excluding header)
  - region root inferred as the parent of ``combined/`` from the GFA path
  - ``{region}/{chunk_id}/shasta/ShastaRun_{suffix}`` where suffix is parsed
    from ``combined_merged_assembly_{suffix}.gfa``

Chunks on chromosomes listed in SKIP_CHROMOSOMES (from ``--config`` YAML) are
excluded from anchor loading.

hap1: longest bp path through each assembled chain.
hap2: graph walk preferring anchor-defined siblings of hap1 contigs; haploid
      contigs (no siblings) appear in both haplotypes.
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# Reuse GFA / chain pipeline from pseudohaps.py
_SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPT_DIR))
import pseudohaps as ph  # noqa: E402

# Reuse Shasta anchor parsing from anchor_stitch.py
sys.path.insert(0, str(_SCRIPT_DIR.parent / "tailor" / "scripts"))
import anchor_stitch  # noqa: E402


# ---------------------------------------------------------------------------
# GFA node <-> Shasta contig naming
# ---------------------------------------------------------------------------

def parse_gfa_node(node_id: str) -> Tuple[Optional[str], str]:
    """'8#0-6-0-0-P0' -> ('8', '0-6-0-0-P0')."""
    if "#" not in node_id:
        return None, node_id
    chunk, chain = node_id.split("#", 1)
    return chunk, chain


def gfa_node_name(chunk_id: str, chain: str) -> str:
    return f"{chunk_id}#{chain}"


# ---------------------------------------------------------------------------
# Snarl / anchor index
# ---------------------------------------------------------------------------

@dataclass
class ContigSnarlProfile:
    """Snarl coverage for one GFA contig node."""
    gfa_node: str
    snarls: Set[int] = field(default_factory=set)
    snarl_to_anchor_ids: Dict[int, Set[str]] = field(default_factory=dict)
    anchor_ids: Set[str] = field(default_factory=set)


@dataclass
class SnarlPairStats:
    shared_snarls: List[int]
    identical_anchor_snarls: List[int]
    allele_snarl_siblings: List[int]
    primary_not_shared: List[int]
    sibling_not_shared: List[int]
    fraction_primary: float
    fraction_sibling: float
    fraction_allele_primary: float
    fraction_allele_sibling: float
    match_fraction: float  # anchor_stitch-style (identical + sibling-pair anchors)


@dataclass
class ChunkAnchorData:
    chunk_id: str
    chunk: anchor_stitch.Chunk
    walk_to_snarl: Dict[str, int]
    anchor_id_to_snarl: Dict[str, int]
    profiles: Dict[str, ContigSnarlProfile]  # gfa_node -> profile


def load_walk_to_snarl(extended_tsv: Path) -> Dict[str, int]:
    """Map anchor walk string (AssemblyDetails / AnchorsFromJson) -> snarl_id."""
    walk_to_snarl: Dict[str, int] = {}
    with extended_tsv.open() as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            try:
                snarl_id = int(parts[1])
            except ValueError:
                continue
            anchor_path = parts[4].strip()
            if anchor_path:
                walk_to_snarl[anchor_path] = snarl_id
    return walk_to_snarl


def build_anchor_id_to_snarl(
    anchors_by_id: Dict[str, anchor_stitch.Anchor],
    walk_to_snarl: Dict[str, int],
) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for aid, anc in anchors_by_id.items():
        sid = walk_to_snarl.get(anc.walk)
        if sid is not None:
            out[aid] = sid
    return out


def build_contig_snarl_profile(
    gfa_node: str,
    chunk_id: str,
    chain: str,
    chunk: anchor_stitch.Chunk,
    anchor_id_to_snarl: Dict[str, int],
) -> Optional[ContigSnarlProfile]:
    contig = chunk.contigs.get(chain)
    if contig is None:
        return None
    profile = ContigSnarlProfile(gfa_node=gfa_node)
    for ca in contig.anchors:
        profile.anchor_ids.add(ca.anchor_id)
        sid = anchor_id_to_snarl.get(ca.anchor_id)
        if sid is None:
            continue
        profile.snarls.add(sid)
        profile.snarl_to_anchor_ids.setdefault(sid, set()).add(ca.anchor_id)
    return profile


def load_chunk_anchor_data(
    chunk_id: str,
    shasta_run_dir: Path,
    extended_tsv: Optional[Path] = None,
) -> ChunkAnchorData:
    if extended_tsv is None:
        # Default: ../anchors/subgraph.anchors.json.subgraph.sizes.extended.tsv
        chunk_root = shasta_run_dir.parent.parent
        extended_tsv = chunk_root / "anchors" / "subgraph.anchors.json.subgraph.sizes.extended.tsv"
    if not extended_tsv.exists():
        raise FileNotFoundError(f"Extended snarl TSV not found: {extended_tsv}")

    chunk = anchor_stitch.load_chunk(chunk_id, shasta_run_dir)
    walk_to_snarl = load_walk_to_snarl(extended_tsv)
    anchor_id_to_snarl = build_anchor_id_to_snarl(chunk.anchors_by_id, walk_to_snarl)

    profiles: Dict[str, ContigSnarlProfile] = {}
    for chain in chunk.contigs:
        gfa_node = gfa_node_name(chunk_id, chain)
        prof = build_contig_snarl_profile(
            gfa_node, chunk_id, chain, chunk, anchor_id_to_snarl
        )
        if prof is not None:
            profiles[gfa_node] = prof

    return ChunkAnchorData(
        chunk_id=chunk_id,
        chunk=chunk,
        walk_to_snarl=walk_to_snarl,
        anchor_id_to_snarl=anchor_id_to_snarl,
        profiles=profiles,
    )


def build_fusion_profile(
    fusion_node: str,
    chunk_data: Dict[str, ChunkAnchorData],
) -> Optional[ContigSnarlProfile]:
    """Build a snarl profile for a fusion contig (name containing ``+``) by
    unioning the profiles of its component contigs.

    A fusion node such as ``0#0-11-3-0-P2+1#0-2-0-1-P2`` is the stitched join of
    one contig per chunk: ``0#0-11-3-0-P2`` (chunk 0) and ``1#0-2-0-1-P2``
    (chunk 1). Each component already has a per-chunk profile, so we combine
    them and let the fusion participate in sibling detection like any other
    contig.

    Snarl ids are namespaced by chunk (``"<chunk>:<snarl>"``) because snarl
    numbering is per-chunk: without namespacing, a snarl id from chunk 0 and an
    unrelated one from chunk 1 could collide. Two fusion alleles built the same
    way are therefore still compared on a consistent, collision-free snarl space.
    """
    profile = ContigSnarlProfile(gfa_node=fusion_node)
    found = False
    for component in fusion_node.split("+"):
        cid, _ = parse_gfa_node(component)
        cdata = chunk_data.get(cid)
        if cdata is None:
            continue
        comp_prof = cdata.profiles.get(component)
        if comp_prof is None:
            continue
        found = True
        profile.anchor_ids |= comp_prof.anchor_ids
        for sid in comp_prof.snarls:
            key = f"{cid}:{sid}"
            profile.snarls.add(key)
            profile.snarl_to_anchor_ids.setdefault(key, set()).update(
                comp_prof.snarl_to_anchor_ids.get(sid, set())
            )
    return profile if found else None


def attach_fusion_profiles(
    chunk_data: Dict[str, ChunkAnchorData],
    gfa_nodes: Set[str],
) -> List[str]:
    """Profile every fusion contig in ``gfa_nodes`` (names containing ``+``) and
    attach it to its first chunk's profile map, so :func:`build_sibling_index`
    compares it against the other contigs that branch from the same junction.

    Returns the list of fusion nodes that were successfully profiled. A fusion
    is grouped under its first chunk; its sibling fusion shares that same
    chunk-boundary prefix, so the pair lands in the same comparison group.
    """
    added: List[str] = []
    for node in gfa_nodes:
        if "+" not in node:
            continue
        cid0, _ = parse_gfa_node(node)
        cdata0 = chunk_data.get(cid0)
        if cdata0 is None or node in cdata0.profiles:
            continue
        prof = build_fusion_profile(node, chunk_data)
        if prof is not None:
            cdata0.profiles[node] = prof
            added.append(node)
    return added


def compare_snarl_pair(
    primary: ContigSnarlProfile,
    sibling: ContigSnarlProfile,
    chunk: anchor_stitch.Chunk,
) -> SnarlPairStats:
    """Compare two contigs at snarl resolution."""
    shared: List[int] = []
    identical: List[int] = []
    allele: List[int] = []
    for sid in sorted(primary.snarls):
        if sid not in sibling.snarls:
            continue
        p_aids = primary.snarl_to_anchor_ids.get(sid, set())
        s_aids = sibling.snarl_to_anchor_ids.get(sid, set())
        if p_aids & s_aids:
            shared.append(sid)
            identical.append(sid)
        elif p_aids and s_aids:
            # Different anchors in the same snarl -> true diploid bubble alleles
            shared.append(sid)
            allele.append(sid)

    primary_not = sorted(primary.snarls - set(shared))
    sibling_not = sorted(sibling.snarls - set(shared))
    n_p = len(primary.snarls)
    n_s = len(sibling.snarls)
    frac_p = len(shared) / n_p if n_p else 0.0
    frac_s = len(shared) / n_s if n_s else 0.0
    frac_allele_p = len(allele) / n_p if n_p else 0.0
    frac_allele_s = len(allele) / n_s if n_s else 0.0

    _, chain_p = parse_gfa_node(primary.gfa_node)
    _, chain_s = parse_gfa_node(sibling.gfa_node)
    # The anchor_stitch match fraction is a single-chunk metric. It does not
    # apply to fusion contigs (which span chunks; their parsed "chain" is not a
    # chunk-local contig) or to any chain missing from this chunk, so fall back
    # to the snarl-fraction criterion (match_fraction = 0.0) in those cases.
    if chain_p in chunk.contigs and chain_s in chunk.contigs:
        match_fraction = anchor_stitch.compare_contigs_for_sibling(
            chunk, chain_p, chain_s
        ).match_fraction
    else:
        match_fraction = 0.0

    return SnarlPairStats(
        shared_snarls=shared,
        identical_anchor_snarls=identical,
        allele_snarl_siblings=allele,
        primary_not_shared=primary_not,
        sibling_not_shared=sibling_not,
        fraction_primary=frac_p,
        fraction_sibling=frac_s,
        fraction_allele_primary=frac_allele_p,
        fraction_allele_sibling=frac_allele_s,
        match_fraction=match_fraction,
    )


def are_sibling_contigs(
    stats: SnarlPairStats,
    *,
    min_match_fraction: float,
    min_shared_snarls: int,
    min_snarl_fraction: float,
) -> bool:
    """
    Two contigs are siblings if they share enough snarl space.

    Key rule: use the *smaller* contig's shared-snarl fraction. A large spine
    contig spans many small alternate-path contigs, so its own shared fraction
    is necessarily low (it has lots of other snarls matching the other small
    contigs). Each small contig, however, sits almost entirely inside the spine.
    The smaller contig (fewer snarls) always has the higher shared fraction, so
    we take max(fraction_primary, fraction_sibling).
    """
    if len(stats.shared_snarls) < min_shared_snarls:
        return False

    if stats.match_fraction >= min_match_fraction:
        return True

    # Smaller contig's shared-snarl fraction (= the higher of the two fractions).
    smaller_frac = max(stats.fraction_primary, stats.fraction_sibling)
    if smaller_frac >= min_snarl_fraction:
        return True

    return False


def build_sibling_index(
    chunk_data: Dict[str, ChunkAnchorData],
    gfa_nodes: Set[str],
    *,
    min_match_fraction: float = 0.5,
    min_shared_snarls: int = 2,
    min_snarl_fraction: float = 0.5,
) -> Tuple[Dict[str, List[str]], List[Tuple[str, str, SnarlPairStats]]]:
    """
    Returns (sibling_dict, pair_records) for all sibling pairs among gfa_nodes.
    sibling_dict maps each contig -> sorted list of its siblings (bidirectional).
    """
    sibling_dict: Dict[str, List[str]] = defaultdict(list)
    pair_records: List[Tuple[str, str, SnarlPairStats]] = []

    # Group nodes by chunk (siblings only within chunk)
    by_chunk: Dict[str, List[str]] = defaultdict(list)
    for node in sorted(gfa_nodes):
        cid, _ = parse_gfa_node(node)
        if cid is not None:
            by_chunk[cid].append(node)

    for cid, nodes in by_chunk.items():
        cdata = chunk_data.get(cid)
        if cdata is None:
            continue
        profiles = cdata.profiles

        for i, a in enumerate(nodes):
            if a not in profiles:
                continue
            for b in nodes[i + 1 :]:
                if b not in profiles:
                    continue
                stats = compare_snarl_pair(profiles[a], profiles[b], cdata.chunk)
                if not are_sibling_contigs(
                    stats,
                    min_match_fraction=min_match_fraction,
                    min_shared_snarls=min_shared_snarls,
                    min_snarl_fraction=min_snarl_fraction,
                ):
                    continue
                pair_records.append((a, b, stats))
                sibling_dict[a].append(b)
                sibling_dict[b].append(a)

    for k in sibling_dict:
        sibling_dict[k] = sorted(set(sibling_dict[k]))
    return dict(sibling_dict), pair_records


def load_sibling_dict_from_tsv(path: str) -> Dict[str, List[str]]:
    """Read a ``*_sibling_contigs.tsv`` (written by :func:`write_sibling_pairs_tsv`)
    back into a bidirectional sibling map ``contig -> sorted siblings``.

    Only the first two columns (primary_contig, sibling_contig) are used; the
    pairing is symmetrised so either column order yields the same map.
    """
    sib: Dict[str, Set[str]] = defaultdict(set)
    with open(path) as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            a, b = parts[0].strip(), parts[1].strip()
            if a and b:
                sib[a].add(b)
                sib[b].add(a)
    return {k: sorted(v) for k, v in sib.items()}


def write_sibling_pairs_tsv(
    path: str,
    pair_records: List[Tuple[str, str, SnarlPairStats]],
) -> None:
    with open(path, "w") as out:
        out.write(
            "primary_contig\tsibling_contig\tshared_snarls\t"
            "identical_anchor_snarls\tallele_snarl_siblings\t"
            "primary_not_shared_snarls\tsibling_not_shared_snarls\t"
            "fraction_shared_snarls_primary\tfraction_shared_snarls_sibling\t"
            "fraction_allele_snarls_primary\tfraction_allele_snarls_sibling\t"
            "match_fraction\n"
        )
        for primary, sibling, stats in pair_records:
            out.write(
                f"{primary}\t{sibling}\t"
                f"{','.join(str(s) for s in stats.shared_snarls)}\t"
                f"{','.join(str(s) for s in stats.identical_anchor_snarls)}\t"
                f"{','.join(str(s) for s in stats.allele_snarl_siblings)}\t"
                f"{','.join(str(s) for s in stats.primary_not_shared)}\t"
                f"{','.join(str(s) for s in stats.sibling_not_shared)}\t"
                f"{stats.fraction_primary:.4f}\t{stats.fraction_sibling:.4f}\t"
                f"{stats.fraction_allele_primary:.4f}\t{stats.fraction_allele_sibling:.4f}\t"
                f"{stats.match_fraction:.4f}\n"
            )
        # Also emit reverse orientation for convenience (primary <-> sibling)
        for primary, sibling, stats in pair_records:
            rev = SnarlPairStats(
                shared_snarls=stats.shared_snarls,
                identical_anchor_snarls=stats.identical_anchor_snarls,
                allele_snarl_siblings=stats.allele_snarl_siblings,
                primary_not_shared=stats.sibling_not_shared,
                sibling_not_shared=stats.primary_not_shared,
                fraction_primary=stats.fraction_sibling,
                fraction_sibling=stats.fraction_primary,
                fraction_allele_primary=stats.fraction_allele_sibling,
                fraction_allele_sibling=stats.fraction_allele_primary,
                match_fraction=stats.match_fraction,
            )
            out.write(
                f"{sibling}\t{primary}\t"
                f"{','.join(str(s) for s in rev.shared_snarls)}\t"
                f"{','.join(str(s) for s in rev.identical_anchor_snarls)}\t"
                f"{','.join(str(s) for s in rev.allele_snarl_siblings)}\t"
                f"{','.join(str(s) for s in rev.primary_not_shared)}\t"
                f"{','.join(str(s) for s in rev.sibling_not_shared)}\t"
                f"{rev.fraction_primary:.4f}\t{rev.fraction_sibling:.4f}\t"
                f"{rev.fraction_allele_primary:.4f}\t{rev.fraction_allele_sibling:.4f}\t"
                f"{rev.match_fraction:.4f}\n"
            )


# ---------------------------------------------------------------------------
# Pseudohaplotype construction
# ---------------------------------------------------------------------------

def path_bp(graph: ph.GFAGraph, path: List[str]) -> int:
    return sum(len(graph.nodes.get(n, "")) for n in path)


def longest_bp_path(
    graph: ph.GFAGraph,
    paths: List[Tuple[Tuple[str, str], List[str]]],
) -> Optional[List[str]]:
    best: Optional[List[str]] = None
    best_bp = -1
    for _, path in paths:
        bp = path_bp(graph, path)
        if bp > best_bp:
            best_bp = bp
            best = path
    return best


def pick_adjacent_sibling(
    contig: str,
    sibling_dict: Dict[str, List[str]],
    exclude: Set[str],
    reachable: Set[str],
) -> Optional[str]:
    """Pick an anchor sibling of `contig` that is a direct successor of the
    current node (i.e. lies in ``reachable``).

    Requiring graph adjacency is essential: a contig may be an anchor sibling of
    `contig` because it covers the same region (e.g. an alternate that spans a
    whole sub-path), but swapping to it mid-walk would create a chimeric,
    graph-invalid jump. Only a sibling that genuinely branches in parallel from
    the current node (a bubble allele, or a fork/stray off the same junction) is
    a valid allele substitution.
    """
    for s in sibling_dict.get(contig, []):
        if s in exclude:
            continue
        if s in reachable:
            return s
    return None


def reconnects_to_chain(
    start: str,
    chain_nodes: Set[str],
    chain_chunks: Set[Optional[str]],
    graph: ph.GFAGraph,
    visited: Set[str],
    max_depth: int = 20,
) -> bool:
    """True if `start` can reach an unvisited chain node within ``max_depth``
    graph hops, traversing only nodes that belong to the chain's own chunk(s).

    This distinguishes a *local* alternate allele — one that detours off the
    chain and rejoins it a few nodes later (an inner bubble, or a convergent
    source/stray contig that links back into the shared flank) — from an edge
    that simply leaves the chain for an unrelated chromosome or chain. The chunk
    restriction keeps the bounded search from wandering across chromosomes.
    """
    if start in chain_nodes and start not in visited:
        return True
    stack: List[Tuple[str, int]] = [(start, 0)]
    seen: Set[str] = {start}
    while stack:
        node, depth = stack.pop()
        if depth >= max_depth:
            continue
        for to, _, _, _ in graph.out_edges.get(node, []):
            if to in seen:
                continue
            if to in chain_nodes and to not in visited:
                return True
            cid, _ = parse_gfa_node(to)
            if cid not in chain_chunks:
                continue
            seen.add(to)
            stack.append((to, depth + 1))
    return False


def pick_stray_sibling(
    contig: str,
    sibling_dict: Dict[str, List[str]],
    exclude: Set[str],
    hap1_nodes: Set[str],
    graph: ph.GFAGraph,
    chain_nodes: Set[str],
    chain_chunks: Set[Optional[str]],
) -> Optional[str]:
    """Pick an anchor sibling of `contig` that is *not* a graph-adjacent bubble
    allele (so ``pick_adjacent_sibling`` missed it), is not itself on the hap1
    backbone, and rejoins this chain.

    These are alternate alleles the assembler emitted as a separate contig that
    does not branch in parallel from the current node — often a graph source/
    sink or a stub linked to the shared flank on only one side (a convergent
    bubble), e.g. ``22#0-3-0-0-P0`` sitting parallel to the long spine
    ``22#0-5-4-0-P1`` and rejoining only downstream at ``22#0-4-0-0-P1``, or
    ``0#0-4-0-0-P0`` parallel to ``0#0-5-0-0-P0`` and rejoining at
    ``0#0-0-0-0-P0``.

    The anchor-sibling relationship is computed up front from shared anchor/
    snarl space, independent of local graph topology, and is the authority for
    what counts as an allele. We still require the candidate to reconnect to the
    chain (directly or within a short, same-chunk detour) so an unrelated
    far-away sibling is not spliced in; once emitted, the caller follows it back
    to the backbone (or resumes hap1 if it dead-ends).
    """
    for s in sibling_dict.get(contig, []):
        if s in exclude or s in hap1_nodes:
            continue
        if reconnects_to_chain(s, chain_nodes, chain_chunks, graph, exclude):
            return s
    return None


def construct_hap2_path(
    graph: ph.GFAGraph,
    hap1: List[str],
    chain_nodes: Set[str],
    sibling_dict: Dict[str, List[str]],
    end_nodes: Set[str],
) -> List[str]:
    """
    Build hap2 by walking the graph while substituting anchor siblings for the
    alleles taken by hap1. hap1 is the backbone guide; at every step the walk is
    in one of two modes:

    * **On the backbone** (``current`` is a hap1 node). Whenever the next hap1
      node has an anchor sibling that hap2 has not used, swap to that sibling so
      the two haplotypes diverge at the bubble. The sibling may be:
        - a *divergent* bubble allele branching in parallel from the current
          node (e.g. ``...-3-1-P2`` vs ``...-3-0-P2``), or
        - a *stray/convergent* allele the assembler split into its own contig
          (often a graph source) that rejoins downstream (e.g. ``22#0-3-0-0-P0``
          parallel to ``22#0-5-4-0-P1``; ``0#0-4-0-0-P0`` parallel to
          ``0#0-5-0-0-P0``).

    * **Off the backbone** (``current`` is an alternate allele just emitted).
      The alternate may itself be a multi-node sub-path with its own inner
      bubble (e.g. ``0#0-12-0-0-P1 -> 0#0-12-1-*-P2 -> 0#0-5-0-0-P0``), so we
      follow graph edges that lead back to the chain, applying the same sibling
      swap when we step onto a hap1 node. If nothing reconnects (a true stray
      dead-end), we resume the hap1 backbone just past the locus the alternate
      replaced instead of truncating hap2.

    The anchor-sibling index is the authority for what is an allele; the chunk-
    bounded reconnection check keeps off-backbone detours from wandering into an
    unrelated chain.
    """
    if not hap1:
        return []

    hap1_nodes = set(hap1)
    hap1_index: Dict[str, int] = {}
    for i, nd in enumerate(hap1):
        hap1_index.setdefault(nd, i)

    chain_chunks: Set[Optional[str]] = {parse_gfa_node(n)[0] for n in chain_nodes}

    # Backbone index that an emitted alternate node substitutes for, so we know
    # where to resume hap1 if that alternate dead-ends. Carried forward across a
    # multi-node alternate sub-path.
    alt_backbone: Dict[str, int] = {}

    # Start hap2 on the alternate allele when the chain begins on a bubble: if
    # hap1's first node has an unused, chain-reconnecting sibling, hap2 should
    # open with that sibling (the two haplotypes must diverge at this bubble too;
    # otherwise both start on the same allele and the alternate is dropped).
    start = hap1[0]
    start_alt = pick_stray_sibling(
        start, sibling_dict, set(), hap1_nodes, graph, chain_nodes, chain_chunks
    )
    if start_alt is not None:
        alt_backbone[start_alt] = hap1_index.get(start, 0)
        current = start_alt
    else:
        current = start
    visited: Set[str] = {current}
    path: List[str] = [current]

    safety = len(chain_nodes) * 4 + 20

    def resume_backbone(node: str) -> Optional[str]:
        """Next unvisited hap1 node past the locus `node` stood in for."""
        bpos = alt_backbone.get(node, hap1_index.get(node))
        if bpos is None:
            return None
        j = bpos + 1
        while j < len(hap1):
            if hap1[j] not in visited:
                return hap1[j]
            j += 1
        return None

    def sibling_substitute(to: str, reachable: Set[str]) -> Optional[str]:
        """The anchor sibling hap2 should take instead of hap1 node `to`."""
        alt = pick_adjacent_sibling(to, sibling_dict, visited, reachable)
        if alt is not None:
            return alt
        return pick_stray_sibling(
            to, sibling_dict, visited, hap1_nodes, graph, chain_nodes, chain_chunks
        )

    for _ in range(safety):
        if current in end_nodes and len(path) > 1:
            break

        on_backbone = current in hap1_nodes

        # Direct successors of the current node (graph-adjacent only).
        reachable: Set[str] = {
            to for to, _, _, _ in graph.out_edges.get(current, [])
        }

        # Candidate next nodes. (node, backbone_index, is_alt) where is_alt marks
        # an alternate allele (off hap1) rather than a plain backbone step.
        candidates: List[Tuple[str, Optional[int], bool]] = []
        for to in reachable:
            if to in visited:
                continue
            on_chain = to in chain_nodes
            if on_backbone:
                # Only leave the chain via an anchor-sibling swap (below).
                if not on_chain:
                    continue
            else:
                # Mid-detour: only follow successors that rejoin this chain.
                if not on_chain and not reconnects_to_chain(
                    to, chain_nodes, chain_chunks, graph, visited
                ):
                    continue

            if to in hap1_nodes:
                alt = sibling_substitute(to, reachable)
                if alt is not None:
                    candidates.append((alt, hap1_index.get(to), True))
                    continue
                candidates.append((to, hap1_index.get(to), False))
            else:
                candidates.append((to, hap1_index.get(to), True))

        if candidates:
            # Prefer alternates (the diploid divergence), then longer contigs.
            candidates.sort(
                key=lambda c: (
                    0 if c[2] else 1,
                    -len(graph.nodes.get(c[0], "")),
                    c[0],
                )
            )
            nxt, bidx, is_alt = candidates[0]
            if is_alt and nxt not in hap1_nodes:
                # Remember which backbone locus this alternate stands in for so a
                # later dead-end can resume hap1; carry it across the sub-path.
                if bidx is not None:
                    alt_backbone[nxt] = bidx
                elif current in alt_backbone:
                    alt_backbone[nxt] = alt_backbone[current]
        else:
            # Dead-end (e.g. a stray contig that never reconnects). Resume the
            # hap1 backbone just past the locus this node stood in for.
            nxt = resume_backbone(current)
            if nxt is None:
                break
            # The resume lands on a raw hap1 node. If that node is a bubble
            # allele with an unused sibling, take the alternate instead so hap2
            # still diverges from hap1 here (otherwise both haps keep the same
            # allele and the alternate is dropped from the assembly entirely).
            alt = pick_stray_sibling(
                nxt, sibling_dict, visited, hap1_nodes, graph, chain_nodes, chain_chunks
            )
            if alt is not None:
                alt_backbone[alt] = hap1_index.get(nxt, len(hap1))
                nxt = alt

        path.append(nxt)
        visited.add(nxt)
        current = nxt

    return path


def construct_pseudohaplotypes_anchors(
    graph: ph.GFAGraph,
    assembled_chains: Dict[str, Dict[int, List[Tuple[Tuple[str, str], List[str]]]]],
    sibling_dict: Dict[str, List[str]],
) -> Dict[str, Dict[int, Tuple[List[str], List[str]]]]:
    """hap1 = longest bp path; hap2 = sibling-aware walk."""
    out: Dict[str, Dict[int, Tuple[List[str], List[str]]]] = defaultdict(dict)

    for chrom, chains in assembled_chains.items():
        for cid, plist in chains.items():
            if not plist:
                continue

            hap1 = longest_bp_path(graph, plist)
            if not hap1:
                print(f"Warning: no hap1 path for {chrom}:{cid}")
                continue

            chain_nodes: Set[str] = set()
            starts: List[str] = []
            ends: List[str] = []
            for (s, e), path in plist:
                chain_nodes.update(path)
                starts.append(s)
                ends.append(e)

            end_nodes = set(ends) | (graph.get_potential_end_nodes() & chain_nodes)
            hap2 = construct_hap2_path(graph, hap1, chain_nodes, sibling_dict, end_nodes)
            if not hap2:
                print(f"Warning: empty hap2 for {chrom}:{cid}; using hap1")
                hap2 = hap1.copy()

            out[chrom][cid] = (hap1, hap2)

    return out


def _has_edge(graph: ph.GFAGraph, a: str, b: str) -> bool:
    """True if the GFA has a direct edge a -> b."""
    return any(to == b for to, _, _, _ in graph.out_edges.get(a, []))


def split_path_into_connected_subchains(
    graph: ph.GFAGraph, path: List[str]
) -> List[List[str]]:
    """Split a haplotype path at fake junctions into graph-connected sub-chains.

    hap2 may step from one node to another that is *not* a real GFA neighbour:
    when a stray/convergent sibling allele is substituted in, or when the walk
    resumes the hap1 backbone after a dead-end. Concatenating across such a
    non-edge would fabricate a junction (and, because per-node overlaps are
    resolved against the real edges only, would emit an incorrect sequence). We
    therefore break the path at every consecutive pair with no edge, so each
    returned sub-chain is a contiguous, graph-valid walk that can be emitted as
    its own FASTA record.

    hap1 is a real start-to-end walk and so always yields a single sub-chain.
    """
    if not path:
        return []
    subchains: List[List[str]] = [[path[0]]]
    for prev, nd in zip(path, path[1:]):
        if _has_edge(graph, prev, nd):
            subchains[-1].append(nd)
        else:
            subchains.append([nd])
    return subchains


def write_pseudohaplotypes_tsv(
    out_path: str,
    pseudohaplotypes: Dict[str, Dict[int, Tuple[List[str], List[str]]]],
    graph: ph.GFAGraph,
) -> None:
    """Write the pseudohaplotypes TSV, one row per graph-connected sub-chain.

    hap1 is always a single sub-chain (subchain_id 1); hap2 may break into
    several (subchain_id 1, 2, ...) at fake junctions.
    """
    with open(out_path, "w") as out:
        out.write(
            "chromosome\tassembled_chain_id\thaplotype_id\tsubchain_id\t"
            "start\tend\tpath_bp\tpath\n"
        )
        for chrom in sorted(pseudohaplotypes):
            for cid in sorted(pseudohaplotypes[chrom]):
                hap1, hap2 = pseudohaplotypes[chrom][cid]
                for hid, full in ((1, hap1), (2, hap2)):
                    subs = split_path_into_connected_subchains(graph, full)
                    for sidx, sub in enumerate(subs, start=1):
                        out.write(
                            f"{chrom}\t{cid}\t{hid}\t{sidx}\t{sub[0]}\t{sub[-1]}\t"
                            f"{path_bp(graph, sub)}\t{','.join(sub)}\n"
                        )


def write_haplotype_fastas(
    output_prefix: str,
    pseudohaplotypes: Dict[str, Dict[int, Tuple[List[str], List[str]]]],
    graph: ph.GFAGraph,
) -> Tuple[str, str]:
    """Write hap1/hap2 FASTAs, one record per graph-connected sub-chain.

    Record header: ``>{chrom}_chain_{cid}_hap{hid}_{subchain}_start_{s}_end_{e}``.
    hap1 emits a single record per chain (subchain 1); hap2 emits one record per
    contiguous sub-chain so no sequence is glued across a non-existent edge.
    """
    hap1_fasta = f"{output_prefix}_hap1.fasta"
    hap2_fasta = f"{output_prefix}_hap2.fasta"
    with open(hap1_fasta, "w") as h1, open(hap2_fasta, "w") as h2:
        for chrom in sorted(pseudohaplotypes):
            for cid in sorted(pseudohaplotypes[chrom]):
                hap1, hap2 = pseudohaplotypes[chrom][cid]
                for hid, full, fh in ((1, hap1, h1), (2, hap2, h2)):
                    subs = split_path_into_connected_subchains(graph, full)
                    for sidx, sub in enumerate(subs, start=1):
                        seq = graph.get_path_sequence(sub)
                        fh.write(
                            f">{chrom}_chain_{cid}_hap{hid}_{sidx}_"
                            f"start_{sub[0]}_end_{sub[-1]}\n"
                        )
                        for i in range(0, len(seq), 80):
                            fh.write(seq[i : i + 80] + "\n")
    return hap1_fasta, hap2_fasta


def write_bandage_csvs(
    output_prefix: str,
    pseudohaplotypes: Dict[str, Dict[int, Tuple[List[str], List[str]]]],
    graph: ph.GFAGraph,
    sibling_dict: Dict[str, List[str]],
) -> Tuple[str, str]:
    """Write Bandage colour CSVs, one per haplotype.

    Each CSV lists *only* the contigs that actually appear in that haplotype's
    paths (hap1 -> red, hap2 -> green), with their ploidy. A contig that is not
    in a given haplotype is omitted from that file, so loading the CSV onto the
    graph never colours a contig the haplotype does not contain. These are
    regenerated alongside the pseudohaplotype TSV/FASTA so they stay in sync.
    """
    ploidy_rows = compute_ploidy_assignments(pseudohaplotypes, graph, sibling_dict)
    contig_to_ploidy: Dict[str, str] = {}
    for _chrom, _cid, node, ploidy, _ctx, _len in ploidy_rows:
        contig_to_ploidy.setdefault(node, ploidy)

    hap1_contigs: Set[str] = set()
    hap2_contigs: Set[str] = set()
    for chains in pseudohaplotypes.values():
        for hap1, hap2 in chains.values():
            hap1_contigs.update(hap1)
            hap2_contigs.update(hap2)

    hap1_file = f"{output_prefix}_hap1_contigs.bandage.csv"
    hap2_file = f"{output_prefix}_hap2_contigs.bandage.csv"
    for fname, contigs, color in (
        (hap1_file, hap1_contigs, "#bf3030"),
        (hap2_file, hap2_contigs, "#30bf30"),
    ):
        with open(fname, "w") as out:
            out.write("Contig,Ploidy,Color\n")
            for contig in sorted(contigs):
                out.write(f"{contig},{contig_to_ploidy.get(contig, 'unknown')},{color}\n")
    return hap1_file, hap2_file


def compute_ploidy_assignments(
    pseudohaplotypes: Dict[str, Dict[int, Tuple[List[str], List[str]]]],
    graph: ph.GFAGraph,
    sibling_dict: Dict[str, List[str]],
) -> List[Tuple[str, str, str, str, str, int]]:
    rows = []
    for chrom, chains in sorted(pseudohaplotypes.items()):
        for cid, (hap1, hap2) in sorted(chains.items()):
            h1, h2 = set(hap1), set(hap2)
            all_nodes = h1 | h2
            for node in sorted(all_nodes):
                length = len(graph.nodes.get(node, ""))
                has_siblings = len(sibling_dict.get(node, [])) > 0
                if node in h1 and node in h2:
                    ploidy = "haploid"
                    context = "in_both_haplotypes"
                elif not has_siblings:
                    ploidy = "haploid"
                    context = "no_anchor_siblings"
                else:
                    ploidy = "diploid"
                    context = "hap1_only" if node in h1 else "hap2_only"
                rows.append((chrom, str(cid), node, ploidy, context, length))
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Chunk / Shasta directory auto-discovery (Snakemake-compatible)
# ---------------------------------------------------------------------------

def infer_region_dir(gfa_path: Path, region_dir: Optional[Path] = None) -> Path:
    """
    Infer the region root directory containing per-chunk folders (0/, 1/, ...).

    Default: parent of ``combined/`` when the GFA lives under ``.../combined/*.gfa``.
    """
    if region_dir is not None:
        return region_dir.resolve()
    gfa_path = gfa_path.resolve()
    if gfa_path.parent.name == "combined":
        return gfa_path.parent.parent
    raise ValueError(
        f"Cannot infer region directory from GFA path {gfa_path}. "
        "Pass --region-dir explicitly."
    )


def extract_shasta_run_suffix(gfa_path: Path) -> Optional[str]:
    """
    From ``combined_merged_assembly_p0.25_k16_c10.gfa`` return ``p0.25_k16_c10``.
    """
    m = re.match(r"combined_merged_assembly_(.+)\.gfa$", gfa_path.name)
    return m.group(1) if m else None


def discover_shasta_run_dir(
    region_dir: Path,
    chunk_id: str,
    shasta_suffix: Optional[str] = None,
) -> Path:
    """Find ``{region_dir}/{chunk_id}/shasta/ShastaRun_*``."""
    shasta_parent = region_dir / chunk_id / "shasta"
    if not shasta_parent.is_dir():
        raise FileNotFoundError(f"Missing shasta directory: {shasta_parent}")

    candidates = sorted(p for p in shasta_parent.glob("ShastaRun_*") if p.is_dir())
    if not candidates:
        raise FileNotFoundError(f"No ShastaRun_* directory under {shasta_parent}")

    if shasta_suffix is not None:
        preferred = shasta_parent / f"ShastaRun_{shasta_suffix}"
        if preferred.is_dir():
            return preferred
        # Fall back if only one candidate exists
        if len(candidates) == 1:
            return candidates[0]
        names = [p.name for p in candidates]
        raise FileNotFoundError(
            f"Expected {preferred.name} under {shasta_parent}, found: {names}"
        )

    if len(candidates) == 1:
        return candidates[0]
    names = [p.name for p in candidates]
    raise FileNotFoundError(
        f"Multiple ShastaRun_* dirs under {shasta_parent}: {names}. "
        "Use a combined_merged_assembly_<suffix>.gfa input or pass --chunk explicitly."
    )


def discover_extended_tsv(region_dir: Path, chunk_id: str) -> Path:
    path = (
        region_dir / chunk_id / "anchors"
        / "subgraph.anchors.json.subgraph.sizes.extended.tsv"
    )
    if not path.exists():
        raise FileNotFoundError(f"Extended snarl TSV not found: {path}")
    return path


def discover_chunks_from_region(
    region_dir: Path,
    chunk_ids: List[str],
    shasta_suffix: Optional[str] = None,
) -> List[Tuple[str, Path, Path]]:
    """
    Return list of (chunk_id, shasta_run_dir, extended_tsv) for each chunk id
    listed in the chunk-points file.
    """
    specs: List[Tuple[str, Path, Path]] = []
    missing: List[str] = []
    for cid in chunk_ids:
        try:
            shasta_dir = discover_shasta_run_dir(region_dir, cid, shasta_suffix)
            ext_tsv = discover_extended_tsv(region_dir, cid)
            specs.append((cid, shasta_dir, ext_tsv))
        except FileNotFoundError as exc:
            missing.append(f"chunk {cid}: {exc}")
    if missing:
        raise FileNotFoundError(
            "Failed to discover Shasta/anchor paths for some chunks:\n  "
            + "\n  ".join(missing)
        )
    return specs


def load_skip_chromosomes_from_config(config_path: Path) -> Set[str]:
    """Load SKIP_CHROMOSOMES from a Snakemake YAML config file."""
    try:
        import yaml  # type: ignore
    except ImportError:
        yaml = None

    with config_path.open() as f:
        if yaml is not None:
            cfg = yaml.safe_load(f) or {}
            val = cfg.get("SKIP_CHROMOSOMES") or []
            return {str(x) for x in val}

        # Fallback: line-based parse when PyYAML is unavailable
        skips: Set[str] = set()
        for line in f:
            line = line.strip()
            if not line.startswith("SKIP_CHROMOSOMES:"):
                continue
            _, _, rest = line.partition(":")
            rest = rest.strip()
            if rest.startswith("[") and rest.endswith("]"):
                inner = rest[1:-1]
                for part in inner.split(","):
                    part = part.strip().strip("'\"")
                    if part:
                        skips.add(part)
            break
        return skips


def apply_skip_chromosomes(
    chunk_to_chrom: Dict[str, str],
    skip_chromosomes: Set[str],
) -> Tuple[Dict[str, str], List[Tuple[str, str]]]:
    """
    Remove chunks whose chromosome is in skip_chromosomes.

    Returns (filtered_chunk_to_chrom, skipped_as_list_of_(chunk_id, chrom)).
    """
    if not skip_chromosomes:
        return chunk_to_chrom, []
    skipped: List[Tuple[str, str]] = []
    filtered: Dict[str, str] = {}
    for cid, chrom in chunk_to_chrom.items():
        if chrom in skip_chromosomes:
            skipped.append((cid, chrom))
        else:
            filtered[cid] = chrom
    return filtered, skipped


def parse_chunk_specs(specs: List[str]) -> List[Tuple[str, Path]]:
    out: List[Tuple[str, Path]] = []
    for spec in specs:
        if "=" not in spec:
            raise ValueError(f"--chunk expects ID=PATH, got {spec!r}")
        cid, path_str = spec.split("=", 1)
        out.append((cid.strip(), Path(path_str.strip())))
    return out


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
        help="Optional manual Shasta run dir override (repeat per chunk). "
             "If omitted, chunks are discovered from the chunk-points TSV and --region-dir.",
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
    region_dir = infer_region_dir(
        gfa_path,
        Path(args.region_dir) if args.region_dir else None,
    )
    shasta_suffix = extract_shasta_run_suffix(gfa_path)
    chunk_to_chrom_all = ph.parse_chunk_points_file(args.chunk_points)

    skip_chromosomes: Set[str] = set(args.skip_chromosome)
    if args.config:
        skip_chromosomes.update(load_skip_chromosomes_from_config(Path(args.config)))

    chunk_to_chrom, skipped_chunks = apply_skip_chromosomes(
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
        by_chrom: Dict[str, List[str]] = defaultdict(list)
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
    chunk_data: Dict[str, ChunkAnchorData] = {}

    if args.chunk:
        chunk_specs = parse_chunk_specs(args.chunk)
        for cid, shasta_dir in chunk_specs:
            ext = extended_overrides.get(cid) or discover_extended_tsv(region_dir, cid)
            chunk_data[cid] = load_chunk_anchor_data(cid, shasta_dir, ext)
            n_prof = len(chunk_data[cid].profiles)
            print(f"  chunk {cid}: {n_prof} contig profiles from {shasta_dir}")
    else:
        discovered = discover_chunks_from_region(region_dir, chunk_ids, shasta_suffix)
        for cid, shasta_dir, ext_tsv in discovered:
            ext = extended_overrides.get(cid, ext_tsv)
            chunk_data[cid] = load_chunk_anchor_data(cid, shasta_dir, ext)
            n_prof = len(chunk_data[cid].profiles)
            print(f"  chunk {cid}: {n_prof} contig profiles from {shasta_dir}")

    print(f"  Done in {time.time() - t0:.2f}s")

    print("\nLoading GFA...")
    t0 = time.time()
    graph = ph.load_gfa(args.gfa)
    print(f"  {len(graph.all_nodes)} nodes in {time.time() - t0:.2f}s")

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

    # Profile fusion contigs (stitched cross-chunk joins, names with '+') from
    # their component contigs so they participate in sibling detection too.
    fusion_added = attach_fusion_profiles(chunk_data, graph.all_nodes)
    print(f"  profiled {len(fusion_added)} fusion contigs from their components")

    # All GFA nodes that have anchor profiles (for sibling TSV)
    anchor_gfa_nodes = set()
    for cdata in chunk_data.values():
        anchor_gfa_nodes.update(cdata.profiles.keys())
    anchor_gfa_nodes &= graph.all_nodes

    print("\nBuilding anchor-based sibling index...")
    t0 = time.time()
    sibling_dict, pair_records = build_sibling_index(
        chunk_data,
        anchor_gfa_nodes,
        min_match_fraction=args.min_match_fraction,
        min_shared_snarls=args.min_shared_snarls,
        min_snarl_fraction=args.min_snarl_fraction,
    )
    sibling_tsv = f"{args.output_prefix}_sibling_contigs.tsv"
    write_sibling_pairs_tsv(sibling_tsv, pair_records)
    print(f"  {len(pair_records)} undirected sibling pairs")
    print(f"  {sum(1 for v in sibling_dict.values() if v)} contigs with >=1 sibling")
    print(f"  Wrote {sibling_tsv} in {time.time() - t0:.2f}s")

    print("\nConstructing pseudohaplotypes (hap1=longest bp, hap2=sibling walk)...")
    t0 = time.time()
    pseudohaplotypes = construct_pseudohaplotypes_anchors(
        graph, assembled_chains, sibling_dict
    )
    n_pairs = sum(len(v) for v in pseudohaplotypes.values())
    print(f"  {n_pairs} haplotype pairs in {time.time() - t0:.2f}s")

    print("\nWriting pseudohaplotypes TSV...")
    write_pseudohaplotypes_tsv(
        f"{args.output_prefix}_pseudohaplotypes.tsv", pseudohaplotypes, graph
    )

    print("\nWriting FASTA files...")
    hap1_fasta, hap2_fasta = write_haplotype_fastas(
        args.output_prefix, pseudohaplotypes, graph
    )
    print(f"  {hap1_fasta}, {hap2_fasta}")

    print("\nWriting ploidy assignments...")
    ploidy_rows = compute_ploidy_assignments(pseudohaplotypes, graph, sibling_dict)
    with open(f"{args.output_prefix}_ploidy_assignments.tsv", "w") as out:
        out.write("chromosome\tassembled_chain_id\tcontig_id\tploidy\tcontext\tlength\n")
        for row in ploidy_rows:
            out.write("\t".join(str(x) for x in row) + "\n")

    print("\nWriting Bandage CSVs...")
    b1, b2 = write_bandage_csvs(args.output_prefix, pseudohaplotypes, graph, sibling_dict)
    print(f"  {b1}, {b2}")

    print("\nDone.")


if __name__ == "__main__":
    main()
