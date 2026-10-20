from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .cluster import ClusterResult, complete_linkage_clusters, pairwise_distance_snp_per_mb


@dataclass
class WindowCall:
    chrom: str
    start0: int
    end0: int
    context_start0: int
    context_end0: int
    window_bp: int
    step_bp: int
    maf_min: float
    n_sites: int
    n_sampled_sites: int
    status: str
    cluster: ClusterResult | None

    @property
    def numeric_labels(self):
        return self.cluster.numeric_labels if self.cluster is not None else None


@dataclass
class ScaleBlock:
    chrom: str
    start0: int
    end0: int
    window_bp: int
    maf_min: float
    calls: list[WindowCall]
    n_bridged_windows: int = 0
    gap_intervals: list[tuple[int, int, str]] = field(default_factory=list)

    @property
    def labels(self):
        valid = [call for call in self.calls if call.cluster is not None]
        return max(valid, key=lambda call: call.cluster.assigned_count).numeric_labels if valid else None


@dataclass
class BlockResult:
    block_id: str
    level: str
    parent_id: str
    chrom: str
    start0: int
    end0: int
    labels: list[str]
    numeric_labels: np.ndarray
    cluster: ClusterResult
    n_sites: int
    n_sampled_sites: int
    support_scales: list[int]
    support_mafs: list[float]
    local_pass_fraction: float
    gaps: list[tuple[int, int, str]] = field(default_factory=list)

    @property
    def length_bp(self):
        return self.end0 - self.start0


def _adjusted_rand_score(left, right):
    """Adjusted Rand index without importing scikit-learn during block scans."""
    left = np.asarray(left)
    right = np.asarray(right)
    if left.size != right.size:
        raise ValueError("label vectors must have equal length")
    if left.size < 2:
        return 1.0
    _, left_inverse = np.unique(left, return_inverse=True)
    _, right_inverse = np.unique(right, return_inverse=True)
    n_left = int(left_inverse.max()) + 1
    n_right = int(right_inverse.max()) + 1
    contingency = np.bincount(
        left_inverse * n_right + right_inverse,
        minlength=n_left * n_right,
    ).reshape(n_left, n_right)

    def combinations2(values):
        values = np.asarray(values, dtype=np.float64)
        return float(np.sum(values * (values - 1.0) / 2.0))

    sum_cells = combinations2(contingency)
    sum_left = combinations2(contingency.sum(axis=1))
    sum_right = combinations2(contingency.sum(axis=0))
    total_pairs = left.size * (left.size - 1.0) / 2.0
    expected = sum_left * sum_right / total_pairs
    maximum = 0.5 * (sum_left + sum_right)
    denominator = maximum - expected
    return 1.0 if abs(denominator) < 1e-15 else float((sum_cells - expected) / denominator)


def label_similarity(left, right):
    if left is None or right is None or len(left) != len(right):
        return -1.0
    return _adjusted_rand_score(left, right)


def interval_genotypes(cache, start0, end0, sample_indices, maf_min=0.0, maf_max=0.5, recalc_subset_maf=False):
    lo, hi = cache.bounds(start0, end0)
    if hi <= lo:
        return np.empty((0, len(sample_indices)), dtype=np.int8), np.empty(0, dtype=np.int64), 0
    positions = np.asarray(cache.positions[lo:hi])
    genotype = np.asarray(cache.genotypes[lo:hi][:, sample_indices], dtype=np.int8)
    if recalc_subset_maf:
        called = genotype >= 0
        an = 2 * called.sum(axis=1)
        ac = np.where(called, genotype, 0).sum(axis=1)
        af = np.divide(ac, an, out=np.full(len(an), np.nan), where=an > 0)
        maf = np.minimum(af, 1.0 - af)
    else:
        maf = np.asarray(cache.maf[lo:hi])
    mask = np.isfinite(maf) & (maf + 1e-12 >= maf_min) & (maf <= maf_max + 1e-12)
    return genotype[mask], positions[mask], int(mask.sum())


def call_interval(
    cache,
    start0,
    end0,
    sample_indices,
    sample_names,
    maf_min,
    threshold,
    genotype_model,
    min_sites,
    max_sites,
    min_hap_fraction,
    min_hap_samples,
    maf_max=0.5,
    recalc_subset_maf=False,
):
    genotype, positions, n_sites = interval_genotypes(
        cache, start0, end0, sample_indices, maf_min, maf_max, recalc_subset_maf
    )
    if n_sites < min_sites:
        return None, genotype, positions
    distance = pairwise_distance_snp_per_mb(
        genotype, n_sites, end0 - start0, genotype_model, max_sites, min_compared_sites=min_sites
    )
    cluster = complete_linkage_clusters(
        distance, sample_names, threshold, min_hap_fraction, min_hap_samples
    )
    return cluster, genotype, positions


def scan_windows(
    cache,
    sample_indices,
    sample_names,
    window_bp,
    step_bp,
    maf_min,
    threshold,
    genotype_model,
    min_sites,
    max_sites,
    min_hap_fraction,
    min_hap_samples,
    fine_context_bp=100000,
    fine_max_context_bp=200000,
    region=None,
):
    chrom_length = int(cache.meta["length"])
    region_start, region_end = region if region else (0, chrom_length)
    calls = []
    for start0 in range(region_start, region_end, step_bp):
        end0 = min(region_end, start0 + window_bp)
        if end0 <= start0:
            continue
        context = window_bp if window_bp > 25000 else max(window_bp, fine_context_bp)
        context_start = max(region_start, start0 - max(0, (context - window_bp) // 2))
        context_end = min(region_end, context_start + context)
        context_start = max(region_start, context_end - context)
        cluster, genotype, _ = call_interval(
            cache, context_start, context_end, sample_indices, sample_names, maf_min, threshold,
            genotype_model, min_sites, max_sites, min_hap_fraction, min_hap_samples,
        )
        if cluster is None and window_bp <= 25000 and fine_max_context_bp > context:
            center = (start0 + end0) // 2
            context_start = max(region_start, center - fine_max_context_bp // 2)
            context_end = min(region_end, context_start + fine_max_context_bp)
            context_start = max(region_start, context_end - fine_max_context_bp)
            cluster, genotype, _ = call_interval(
                cache, context_start, context_end, sample_indices, sample_names, maf_min, threshold,
                genotype_model, min_sites, max_sites, min_hap_fraction, min_hap_samples,
            )
        status = "PASS" if cluster is not None and cluster.status == "PASS" else "LOW_INFORMATION" if cluster is None else cluster.status
        calls.append(WindowCall(
            cache.chrom, start0, end0, context_start, context_end, window_bp, step_bp, maf_min,
            int(genotype.shape[0]), min(int(genotype.shape[0]), max_sites), status, cluster,
        ))
        if end0 == region_end:
            break
    return calls


def merge_window_calls(calls, ari_threshold, max_gap_bp, max_gap_fraction, min_block_bp):
    blocks = []
    current = []
    pending = []
    first_labels = last_labels = None
    gaps = []

    def finalize():
        nonlocal current, pending, first_labels, last_labels, gaps
        if current:
            start0 = current[0].start0
            end0 = current[-1].end0
            if end0 - start0 >= min_block_bp:
                blocks.append(ScaleBlock(
                    current[0].chrom, start0, end0, current[0].window_bp, current[0].maf_min,
                    list(current), sum(1 for x in gaps), list(gaps),
                ))
        current, pending, gaps = [], [], []
        first_labels = last_labels = None

    for call in calls:
        if call.status != "PASS":
            if current:
                pending.append(call)
            continue
        if not current:
            current = [call]
            first_labels = last_labels = call.numeric_labels
            pending = []
            continue
        pending_start = pending[0].start0 if pending else current[-1].end0
        pending_end = pending[-1].end0 if pending else pending_start
        # Sliding windows overlap.  Limit the physical interval not covered by
        # the two flanking PASS windows, rather than the full width of every
        # low-information window.
        pending_bp = max(0, call.start0 - current[-1].end0)
        proposed_span = call.end0 - current[0].start0
        fraction_limit = int(max_gap_fraction * proposed_span)
        allowed_gap = min(max_gap_bp, fraction_limit) if max_gap_fraction > 0 else 0
        compatible = (
            label_similarity(last_labels, call.numeric_labels) >= ari_threshold and
            label_similarity(first_labels, call.numeric_labels) >= max(0.0, ari_threshold - 0.10)
        )
        if compatible and pending_bp <= allowed_gap:
            if pending:
                gaps.append((pending_start, pending_end, "BRIDGED_LOW_INFORMATION_OR_DISCORDANT_WINDOWS"))
            current.append(call)
            last_labels = call.numeric_labels
            pending = []
        else:
            finalize()
            current = [call]
            first_labels = last_labels = call.numeric_labels
    finalize()
    return blocks


def merge_scale_blocks(scale_blocks, ari_threshold, max_gap_bp, max_gap_fraction):
    candidates = []
    for block in sorted(scale_blocks, key=lambda b: (b.chrom, b.start0, b.end0, -b.window_bp)):
        if not candidates or candidates[-1]["chrom"] != block.chrom:
            candidates.append({
                "chrom": block.chrom, "start0": block.start0, "end0": block.end0,
                "labels": block.labels, "blocks": [block], "gaps": list(block.gap_intervals),
            })
            continue
        current = candidates[-1]
        gap = max(0, block.start0 - current["end0"])
        proposed = max(current["end0"], block.end0) - current["start0"]
        fraction_limit = int(max_gap_fraction * proposed)
        allowed = min(max_gap_bp, fraction_limit) if max_gap_fraction > 0 else 0
        overlap_or_gap = block.start0 <= current["end0"] + allowed
        compatible = label_similarity(current["labels"], block.labels) >= ari_threshold
        if overlap_or_gap and compatible:
            if gap:
                current["gaps"].append((current["end0"], block.start0, "BRIDGED_BETWEEN_SCALES"))
            current["end0"] = max(current["end0"], block.end0)
            current["blocks"].append(block)
            if block.end0 >= current["end0"]:
                current["labels"] = block.labels
            current["gaps"].extend(block.gap_intervals)
        else:
            candidates.append({
                "chrom": block.chrom, "start0": block.start0, "end0": block.end0,
                "labels": block.labels, "blocks": [block], "gaps": list(block.gap_intervals),
            })
    return candidates


def validate_candidates(
    candidates,
    caches,
    sample_indices,
    sample_names,
    level,
    id_prefix,
    primary_maf,
    threshold,
    genotype_model,
    min_sites,
    max_sites,
    min_hap_fraction,
    min_hap_samples,
    min_scale_support=1,
    parent_id="NA",
):
    results = []
    counter = 0
    for candidate in candidates:
        support_scales = sorted({block.window_bp for block in candidate["blocks"]}, reverse=True)
        support_mafs = sorted({block.maf_min for block in candidate["blocks"]})
        if len(support_scales) < min_scale_support:
            continue
        cache = caches[candidate["chrom"]]
        cluster, genotype, _ = call_interval(
            cache, candidate["start0"], candidate["end0"], sample_indices, sample_names,
            primary_maf, threshold, genotype_model, min_sites, max_sites, min_hap_fraction, min_hap_samples,
        )
        if cluster is None or cluster.status != "PASS":
            continue
        counter += 1
        block_id = f"{candidate['chrom']}_{id_prefix}{counter:04d}"
        local_calls = [call for block in candidate["blocks"] for call in block.calls]
        local_pass = sum(call.status == "PASS" for call in local_calls) / max(1, len(local_calls))
        results.append(BlockResult(
            block_id, level, parent_id, candidate["chrom"], candidate["start0"], candidate["end0"],
            cluster.labels, cluster.numeric_labels, cluster, int(genotype.shape[0]),
            min(int(genotype.shape[0]), max_sites), support_scales, support_mafs, local_pass,
            candidate.get("gaps", []),
        ))
    return results


def fine_children_for_parent(
    parent: BlockResult,
    fine_blocks_by_scale,
    cache,
    sample_indices,
    sample_names,
    primary_maf,
    threshold,
    genotype_model,
    min_sites,
    max_sites,
    min_hap_fraction,
    min_hap_samples,
):
    chosen = []
    for scale in sorted(fine_blocks_by_scale):
        overlapping = [
            block for block in fine_blocks_by_scale[scale]
            if block.chrom == parent.chrom and block.end0 > parent.start0 and block.start0 < parent.end0
        ]
        if overlapping:
            chosen = overlapping
            break
    results = []
    counter = 0
    for block in chosen:
        start0, end0 = max(parent.start0, block.start0), min(parent.end0, block.end0)
        if end0 <= start0:
            continue
        cluster, genotype, _ = call_interval(
            cache, start0, end0, sample_indices, sample_names, primary_maf, threshold,
            genotype_model, min_sites, max_sites, min_hap_fraction, min_hap_samples,
        )
        if cluster is None or cluster.status != "PASS":
            continue
        counter += 1
        results.append(BlockResult(
            f"{parent.block_id}_F{counter:03d}", "FINE", parent.block_id, parent.chrom,
            start0, end0, cluster.labels, cluster.numeric_labels, cluster,
            int(genotype.shape[0]), min(int(genotype.shape[0]), max_sites), [block.window_bp],
            [block.maf_min], sum(call.status == "PASS" for call in block.calls) / max(1, len(block.calls)),
            list(block.gap_intervals),
        ))
    return results
