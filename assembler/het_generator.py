import multiprocessing
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict

from assembler.config import settings

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def find_adjacent_snarl_pairs(snarl_list, snarl_to_anchors_dictionary, read_to_snarl_dictionary):
    """
    This function uses the read journeys to find adjacent snarl pairs.
    """
    adjacent_snarl_pairs = []

    for current_snarl_id in snarl_list:
        current_snarl_anchors = snarl_to_anchors_dictionary[current_snarl_id]
        # print(f"current snarl ID: {current_snarl_id} has {len(current_snarl_anchors)} anchors", flush=True)

        potential_prev_snarl_ids = Counter()
        potential_succ_snarl_ids = Counter()        
        
        for anchor in current_snarl_anchors:
            # print(f" anchor: {anchor!r} has {len(anchor.bp_matched_reads)} reads", flush=True)
            for read_id, strand, *_ in anchor.bp_matched_reads:
                anchor_rank = anchor.read_ranks[read_id]
                read_journey_snarls = read_to_snarl_dictionary[read_id]
                # print(f"  read ID: {read_id}, strand: {strand}, anchor rank: {anchor_rank}", flush=True)
                # print(f"  read journey snarls: {read_journey_snarls}", flush=True)
                
                # snarl before the current anchor
                prev_snarl = (
                    read_journey_snarls[anchor_rank - 1] if anchor_rank > 0 else "NA"
                )

                # snarl after the current anchor
                succ_snarl = (
                    read_journey_snarls[anchor_rank + 1] if anchor_rank + 1 < len(read_journey_snarls) else "NA"
                )

                if strand == 0:
                    potential_prev_snarl_ids[prev_snarl] += 1
                    potential_succ_snarl_ids[succ_snarl] += 1
                    # print(f"  prev snarl: {prev_snarl}", flush=True)
                    # print(f"  succ snarl: {succ_snarl}", flush=True)
                else:
                    potential_prev_snarl_ids[succ_snarl] += 1
                    potential_succ_snarl_ids[prev_snarl] += 1
                    # print(f"  prev snarl: {succ_snarl}", flush=True)
                    # print(f"  succ snarl: {prev_snarl}", flush=True)

        if not potential_prev_snarl_ids and not potential_succ_snarl_ids:
            continue

        # print(f"   potential prev snarls: {potential_prev_snarl_ids}", flush=True)
        # print(f"   potential succ snarls: {potential_succ_snarl_ids}", flush=True)

        # select the most frequent prev and succ snarls
        most_frequent_prev_snarl_id = (
            potential_prev_snarl_ids.most_common(1)[0][0] if potential_prev_snarl_ids else "NA"
        )
        most_frequent_succ_snarl_id = (
            potential_succ_snarl_ids.most_common(1)[0][0] if potential_succ_snarl_ids else "NA"
        )
    
        if most_frequent_prev_snarl_id != "NA":
            adjacent_snarl_pairs.append((most_frequent_prev_snarl_id, current_snarl_id))
        if most_frequent_succ_snarl_id != "NA":
            adjacent_snarl_pairs.append((current_snarl_id, most_frequent_succ_snarl_id))

        # print(f"   adjacent snarl pairs: {adjacent_snarl_pairs}", flush=True)

    return adjacent_snarl_pairs


# ======================================================================================
# Read-based het-anchor detection (POA of common reads' inter-snarl subsequences).
#
# For each adjacent snarl pair (s1, s2): take the reads common to both snarls, slice each
# read's subsequence spanning the INTER-SNARL GAP (excluding the snarl cores, which differ
# across alleles), run abpoa to build an MSA, call clean 2-allele SNP columns and indel
# events, and emit per site TWO byte-identical allele anchors. Ported from the standalone
# prototype scratch_msa_het_finder.py, adapted to run on in-memory anchor objects.
#
# Coordinate convention: anchor read entries store [read_id, strand, start, end] in the
# read's ALIGNMENT orientation. Per read: oriented = seq if strand==0 else revcomp(seq);
# an anchor occupies oriented[start:end]. Emitted positions are in that same frame.
# ======================================================================================

# Field positions inside an anchor read entry (bp_matched_reads); mirrors config.ini.
_READ_ID, _READ_STRAND, _READ_START, _READ_END = 0, 1, 2, 3


def _revcomp(s):
    return s.translate(_COMPLEMENT)[::-1]


def _run_abpoa_msa(abpoa_bin, named_seqs, tmp_fa):
    """named_seqs: list of (name, seq). Returns {name: aligned_row} from abpoa -r1 MSA.
    Rows are parsed by FASTA header (not order), so reordering is safe."""
    with open(tmp_fa, "w") as f:
        for name, seq in named_seqs:
            f.write(f">{name}\n{seq}\n")
    proc = subprocess.run([abpoa_bin, "-r1", tmp_fa], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"abpoa failed (rc={proc.returncode}): {proc.stderr[-500:]}")
    rows = {}
    name = None
    buf = []
    for line in proc.stdout.splitlines():
        if line.startswith(">"):
            if name is not None:
                rows[name] = "".join(buf)
            name = line[1:].split()[0].strip()
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        rows[name] = "".join(buf)
    return rows


def _call_het_columns(rows_by_read, params):
    """Balanced filter over MSA columns. rows_by_read: {read_name: aligned_row}.
    Returns a list of dicts per het column: {col, base_a0, base_a1, reads_a0, reads_a1}."""
    names = list(rows_by_read.keys())
    depth = len(names)
    if depth < params["min_reads"]:
        return []
    L = len(next(iter(rows_by_read.values())))
    min_frac = params["min_allele_frac"]
    min_cnt = params["min_allele_reads"]
    max_other_frac = params["max_other_frac"]

    het = []
    for c in range(L):
        col = [rows_by_read[n][c] for n in names]
        cnt = Counter(col)
        gap = cnt.get("-", 0)
        acgt = [(cnt.get(b, 0), b) for b in "ACGT"]
        acgt.sort(reverse=True)
        (c1, b1), (c2, b2) = acgt[0], acgt[1]
        other = depth - c1 - c2  # 3rd/4th base + gaps + N
        if (c1 >= min_cnt and c2 >= min_cnt
                and c1 >= min_frac * depth and c2 >= min_frac * depth
                and other <= max_other_frac * depth
                and gap <= max_other_frac * depth):
            reads_a0 = [n for n in names if rows_by_read[n][c] == b1]
            reads_a1 = [n for n in names if rows_by_read[n][c] == b2]
            het.append({"col": c, "base_a0": b1, "base_a1": b2,
                        "reads_a0": reads_a0, "reads_a1": reads_a1})
    return het


def _call_indel_events(rows_by_read, params):
    """Find biallelic indel events: maximal runs of MSA columns where reads split cleanly
    into 'has-gap' vs 'has-base', the SAME reads gapped throughout, non-gap side ~monomorphic.
    Returns list of dicts: {i0, i1, base_reads, gap_reads} (indel spans columns i0..i1)."""
    names = list(rows_by_read.keys())
    depth = len(names)
    if depth < params["min_reads"]:
        return []
    L = len(next(iter(rows_by_read.values())))
    min_frac = params["min_allele_frac"]
    min_cnt = params["min_allele_reads"]

    def col_stats(c):
        gapped = frozenset(n for n in names if rows_by_read[n][c] == "-")
        base_cnt = Counter(rows_by_read[n][c] for n in names if rows_by_read[n][c] != "-")
        return gapped, base_cnt

    def is_clean_indel_col(gapped, base_cnt):
        ng = depth - len(gapped)
        if not (len(gapped) >= min_cnt and len(gapped) >= min_frac * depth
                and ng >= min_cnt and ng >= min_frac * depth):
            return False
        if not base_cnt:
            return False
        return base_cnt.most_common(1)[0][1] >= 0.9 * ng  # non-gap ~single base

    events = []
    c = 0
    while c < L:
        gapped, base_cnt = col_stats(c)
        if is_clean_indel_col(gapped, base_cnt):
            i0 = i1 = c
            cc = c + 1
            while cc < L:  # extend while the SAME reads stay gapped
                g2, bc2 = col_stats(cc)
                if is_clean_indel_col(g2, bc2) and g2 == gapped:
                    i1 = cc
                    cc += 1
                else:
                    break
            base_reads = [n for n in names if n not in gapped]
            gap_reads = [n for n in names if n in gapped]
            events.append({"i0": i0, "i1": i1,
                           "base_reads": base_reads, "gap_reads": gap_reads})
            c = cc
        else:
            c += 1
    return events


def _build_allele_anchor(rows_by_read, names, wL, wR, offset0, strand_of, min_reads,
                         read_anchor_intervals):
    """Byte-identical anchor for one allele group over MSA columns [wL, wR].

    The anchor sequence is the gap-COLLAPSED consensus substring; only reads whose collapsed
    window substring equals it are kept (span identical across the anchor's reads). A read is
    dropped when this new anchor's span would overlap any EXISTING anchor already on it
    (shasta2 rejects the read's journey otherwise). Existing spans are widened to >= 1bp so a
    new anchor cannot start at the exact position of a 0bp existing anchor.

    Returns (anchor_seq, entries) with entries=[[name, strand, start, end], ...], or None if
    the anchor would be <2bp or fewer than min_reads survive."""
    if not names:
        return None
    subs = {n: rows_by_read[n][wL:wR + 1].replace("-", "") for n in names}
    cons, _ = Counter(subs.values()).most_common(1)[0]
    if len(cons) < 2:
        return None
    entries = []
    for n in names:
        if subs[n] != cons:               # byte-identity within the allele
            continue
        read_start = offset0[n] + (wL - rows_by_read[n][:wL].count("-"))
        read_end = read_start + len(cons)
        st = strand_of[n]
        if any(lo < read_end and read_start < max(hi, lo + 1)
               for lo, hi in read_anchor_intervals.get((n, st), ())):
            continue                      # would overlap an existing anchor on this read
        entries.append([n, st, read_start, read_end])
    if len(entries) < min_reads:
        return None
    return cons, entries


def _process_snarl_pair(s1, s2, snarl_reads, read_seqs, params, tmp_fa, sites_out, stats,
                        read_anchor_intervals):
    """Find read-based hets in the inter-snarl gap between adjacent snarls s1, s2. Appends
    accepted site dicts to sites_out."""
    dA = snarl_reads.get(s1)
    dB = snarl_reads.get(s2)
    stats["intervals"] += 1
    if not dA or not dB:
        stats["missing_snarl"] += 1
        return

    common = sorted(set(dA) & set(dB))    # sorted -> deterministic order into abpoa
    named_seqs = []
    offset0 = {}       # read_name -> oriented index of first char in its gap slice
    strand_of = {}
    region_lens = []
    for name in common:
        sA, aLo, aHi = dA[name]
        sB, bLo, bHi = dB[name]
        if sA != sB:
            continue
        seq = read_seqs.get(name)
        if seq is None:
            continue
        oriented = seq if sA == 0 else _revcomp(seq)
        # Inter-snarl GAP only (exclude the snarl cores, which differ across alleles).
        if aHi <= bLo:               # s1 left of s2
            lo, hi = aHi, bLo
        elif bHi <= aLo:             # s2 left of s1
            lo, hi = bHi, aLo
        else:
            continue                 # anchor cores overlap on this read
        if hi - lo < params["min_gap"]:
            continue
        named_seqs.append((name, oriented[lo:hi]))
        offset0[name] = lo
        strand_of[name] = sA
        region_lens.append(hi - lo)

    if len(named_seqs) < params["min_reads"]:
        return

    med_len = sorted(region_lens)[len(region_lens) // 2]
    if params["max_interval_bp"] and med_len > params["max_interval_bp"]:
        stats["too_long"] += 1
        return

    try:
        msa = _run_abpoa_msa(params["abpoa_bin"], named_seqs, tmp_fa)
    except RuntimeError as e:
        stats["abpoa_error"] += 1
        print(f"[warn] read-based het: snarl {s1}->{s2}: {e}", file=sys.stderr)
        return

    rows_by_read = {n: msa[n] for n, _ in named_seqs if n in msa}
    msa_len = len(next(iter(rows_by_read.values()))) if rows_by_read else 0
    het_cols = _call_het_columns(rows_by_read, params)
    indel_events = [] if not params["call_indels"] else _call_indel_events(rows_by_read, params)

    # A candidate site = a window [wL,wR] of MSA columns + a 2-group read partition. Each
    # allele anchor is the gap-collapsed byte-identical consensus; a site is kept only if
    # BOTH alleles are valid and distinct. (SNP -> 3bp/3bp; indel -> L+ins+R vs 2bp L+R.)
    candidates = []
    for hc in het_cols:
        c = hc["col"]
        candidates.append({"wL": c - 1, "wR": c + 1,
                           "groups": [hc["reads_a0"], hc["reads_a1"]]})
    for ev in indel_events:
        candidates.append({"wL": ev["i0"] - 1, "wR": ev["i1"] + 1,
                           "groups": [ev["base_reads"], ev["gap_reads"]]})

    built_sites = []
    for cand in candidates:
        wL, wR = cand["wL"], cand["wR"]
        if wL < 0 or wR >= msa_len:
            stats["site_skipped"] += 1
            continue
        anchors, ok = [], True
        for grp_names in cand["groups"]:
            res = _build_allele_anchor(rows_by_read, grp_names, wL, wR, offset0, strand_of,
                                       params["min_anchor_reads"], read_anchor_intervals)
            if res is None:
                ok = False
                break
            anchors.append(res)  # (seq, entries)
        if not ok or anchors[0][0] == anchors[1][0]:
            stats["site_skipped"] += 1
            continue
        built_sites.append({"wL": wL, "wR": wR, "anchors": anchors,
                            "support": sum(len(a[1]) for a in anchors)})

    # Non-overlap resolution within this pair: two anchors on a read must not overlap. Sites
    # whose MSA windows overlap (even sharing one flank column = 1bp read overlap) conflict;
    # greedily keep the higher-support site and drop the rest.
    built_sites.sort(key=lambda s: (-s["support"], s["wL"]))
    accepted = []
    for s in built_sites:
        if any(not (s["wR"] < a["wL"] or a["wR"] < s["wL"]) for a in accepted):
            stats["overlap_dropped"] += 1
            continue
        accepted.append(s)

    if accepted:
        stats["with_variant"] += 1
    accepted.sort(key=lambda s: s["wL"])   # left-to-right, stable
    for s in accepted:
        (seq0, entries0), (seq1, entries1) = s["anchors"]
        variant_type = "SNP" if len(seq0) == len(seq1) else "INDEL"
        sites_out.append({
            "s1": s1, "s2": s2, "variant_type": variant_type,
            "alleles": [entries0, entries1], "seqs": [seq0, seq1],
        })
        stats["sites_emitted"] += 1
        stats["snp_emitted" if variant_type == "SNP" else "indel_emitted"] += 1


def _resolve_site_overlaps_by_coverage(sites):
    """Global overlap resolution across het SITES from different pairs (branching adjacencies
    -> overlapping gaps -> two pairs emit an overlapping het). Two sites whose anchor spans
    overlap on any read cannot coexist (shasta2 orders each read's anchors). Keep the HIGHER-
    coverage site (distinct (read, strand) across both alleles), drop the lower. A site is
    kept/dropped as a unit. Returns (kept_sites, n_dropped)."""
    cov = []
    for site in sites:
        rs = set()
        sp = []
        for entries in site["alleles"]:
            for n, st, a, b in entries:
                rs.add((n, st))
                sp.append((n, st, a, b))
        cov.append((len(rs), sp))
    occupied = defaultdict(list)
    kept_idx = set()
    for i in sorted(range(len(sites)), key=lambda j: -cov[j][0]):   # higher coverage first
        sp = cov[i][1]
        if any(any(lo < b and a < hi for lo, hi in occupied[(n, st)])
               for n, st, a, b in sp):
            continue                   # overlaps a kept (higher-cov) site -> drop
        kept_idx.add(i)
        for n, st, a, b in sp:
            occupied[(n, st)].append((a, b))
    kept = [s for i, s in enumerate(sites) if i in kept_idx]
    return kept, len(sites) - len(kept)


def _new_stats():
    return {"intervals": 0, "with_variant": 0, "sites_emitted": 0, "site_skipped": 0,
            "snp_emitted": 0, "indel_emitted": 0, "overlap_dropped": 0,
            "missing_snarl": 0, "too_long": 0, "abpoa_error": 0}


# --- parallel POA over snarl pairs (embarrassingly parallel; only READS shared data) ------
# The heavy per-pair work (abpoa MSA + het calling) is independent across pairs and touches
# only read-only inputs (snarl_reads, read_sequences, read_anchor_intervals, params). We hand
# those to the workers via a module global set BEFORE forking the pool, so fork's copy-on-write
# gives access without pickling the (large) read_sequences. Each pair's per-read overlap check
# against EXISTING anchors happens here (immutable read_anchor_intervals); the global
# cross-site overlap resolution runs afterwards in the parent — neither is affected by the
# parallelism. Results are gathered in input order (pool.map) so the outcome is deterministic.
_RBH_SHARED = None


def _rbh_worker_init():
    # Nothing to do: on fork the worker already sees _RBH_SHARED via copy-on-write.
    pass


def _process_pair_chunk(chunk_pairs):
    shared = _RBH_SHARED
    params = shared["params"]
    tmp_fa = os.path.join(params["tmp_dir"], f"abpoa_msa_het.{os.getpid()}.fa")
    sites, stats = [], _new_stats()
    for (s1, s2) in chunk_pairs:
        _process_snarl_pair(s1, s2, shared["snarl_reads"], shared["read_sequences"], params,
                            tmp_fa, sites, stats, shared["read_anchor_intervals"])
    try:
        os.remove(tmp_fa)
    except OSError:
        pass
    return sites, stats


def _run_pairs_serial(pairs, snarl_reads, read_sequences, read_anchor_intervals, params):
    tmp_fa = os.path.join(params["tmp_dir"], f"abpoa_msa_het.{os.getpid()}.fa")
    sites, stats = [], _new_stats()
    for (s1, s2) in pairs:
        _process_snarl_pair(s1, s2, snarl_reads, read_sequences, params, tmp_fa, sites, stats,
                            read_anchor_intervals)
    try:
        os.remove(tmp_fa)
    except OSError:
        pass
    return sites, stats


def _run_pairs_parallel(pairs, snarl_reads, read_sequences, read_anchor_intervals, params,
                        threads):
    global _RBH_SHARED
    _RBH_SHARED = {"snarl_reads": snarl_reads, "read_sequences": read_sequences,
                   "read_anchor_intervals": read_anchor_intervals, "params": params}
    # Use more chunks than workers so the pool load-balances (per-pair POA time varies widely).
    n_chunks = min(len(pairs), threads * 4)
    chunk_size = (len(pairs) + n_chunks - 1) // n_chunks
    chunks = [pairs[i:i + chunk_size] for i in range(0, len(pairs), chunk_size)]
    try:
        with multiprocessing.Pool(processes=threads, initializer=_rbh_worker_init) as pool:
            results = pool.map(_process_pair_chunk, chunks, chunksize=1)   # ordered -> deterministic
    finally:
        _RBH_SHARED = None
    sites, stats = [], _new_stats()
    for chunk_sites, chunk_stats in results:
        sites.extend(chunk_sites)
        for k, v in chunk_stats.items():
            stats[k] += v
    return sites, stats


def find_read_based_hets(adjacent_snarl_pairs, snarl_to_anchors_dictionary, read_sequences,
                         params, log=None):
    """Discover read-based het (SNP/indel) sites in the gaps between adjacent snarls.

    Parameters
    ----------
    adjacent_snarl_pairs : list of (s1, s2) tuples (may contain duplicates / both directions)
    snarl_to_anchors_dictionary : {snarl_id: [AnchorBase, ...]} — reads read from bp_matched_reads
    read_sequences : {read_id: sequence_str}
    params : dict with keys abpoa_bin, tmp_dir, min_reads, min_gap, min_allele_frac,
             min_allele_reads, max_other_frac, min_anchor_reads, call_indels, max_interval_bp
    log : optional callable(str) for progress messages

    Returns
    -------
    list of site dicts: {"s1", "s2", "variant_type" ("SNP"/"INDEL"),
                         "alleles": [entries_a0, entries_a1], "seqs": [seq0, seq1]}
    where entries = [[read_id, strand, start, end], ...] in read alignment orientation.
    """
    if log is None:
        def log(_msg):
            return None
    t0 = time.time()

    # 1. Dedup pairs -> sorted unique (min, max) int tuples.
    uniq = set()
    for p in adjacent_snarl_pairs:
        a, b = int(p[0]), int(p[1])
        if a != b:
            uniq.add((min(a, b), max(a, b)))
    pairs = sorted(uniq)

    # 2. snarl_reads: {snarl_id: {read_id: (strand, lo, hi)}} for snarls appearing in pairs,
    #    aggregating reads across ALL allele anchors of the snarl (keep first occurrence).
    needed_snarls = set()
    for a, b in pairs:
        needed_snarls.add(a)
        needed_snarls.add(b)
    snarl_reads = {}
    for sid in needed_snarls:
        d = {}
        for anchor in snarl_to_anchors_dictionary.get(sid, []):
            for read in anchor.bp_matched_reads:
                name = read[_READ_ID]
                if name not in d:
                    d[name] = (read[_READ_STRAND],
                               min(read[_READ_START], read[_READ_END]),
                               max(read[_READ_START], read[_READ_END]))
        if d:
            snarl_reads[sid] = d
    needed_reads = set()
    for d in snarl_reads.values():
        needed_reads.update(d.keys())

    # 3. Existing-anchor spans per (read, strand), over ALL anchors, restricted to the reads
    #    that could carry a new het — so injected het-anchors never overlap an anchor already
    #    on a read (shasta2 would reject that read's journey otherwise).
    read_anchor_intervals = defaultdict(list)
    for anchors in snarl_to_anchors_dictionary.values():
        for anchor in anchors:
            for read in anchor.bp_matched_reads:
                name = read[_READ_ID]
                if name in needed_reads:
                    lo = min(read[_READ_START], read[_READ_END])
                    hi = max(read[_READ_START], read[_READ_END])
                    if lo == hi:      # 0bp anchor: occupy 1bp so a new het can't start on it
                        hi = lo + 1
                    read_anchor_intervals[(name, read[_READ_STRAND])].append((lo, hi))

    log(f"read-based het: {len(pairs)} unique adjacent pairs, "
        f"{len(snarl_reads)} snarls with reads, {len(needed_reads)} reads "
        f"(setup {time.time() - t0:.1f}s)")

    os.makedirs(params["tmp_dir"], exist_ok=True)
    threads = int(params.get("threads", 1) or 1)

    # Per-pair POA is embarrassingly parallel; fan it out across `threads` workers. The global
    # cross-site overlap resolution below then runs once, serially, in the parent.
    t1 = time.time()
    if threads > 1 and len(pairs) > 1:
        sites, stats = _run_pairs_parallel(pairs, snarl_reads, read_sequences,
                                           read_anchor_intervals, params, threads)
    else:
        sites, stats = _run_pairs_serial(pairs, snarl_reads, read_sequences,
                                         read_anchor_intervals, params)

    sites, n_ov = _resolve_site_overlaps_by_coverage(sites)

    log(f"read-based het: intervals={stats['intervals']}, "
        f"intervals_with_variant={stats['with_variant']}, "
        f"sites_emitted={stats['sites_emitted']} "
        f"(snp={stats['snp_emitted']}, indel={stats['indel_emitted']}), "
        f"cross-pair overlap-dropped={n_ov}, kept_sites={len(sites)} "
        f"(POA {time.time() - t1:.1f}s on {threads} thread(s))")
    return sites
