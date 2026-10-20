from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from .common import stable_subsample_indices


@dataclass
class DistanceResult:
    distance: np.ndarray
    compared: np.ndarray
    mismatch: np.ndarray
    n_eligible_sites: int
    n_sampled_sites: int
    span_bp: int


@dataclass
class ClusterResult:
    labels: list[str]
    numeric_labels: np.ndarray
    raw_labels: np.ndarray
    groups: dict[str, list[int]]
    n_major_haps: int
    assigned_count: int
    assigned_fraction: float
    max_within_distance: float
    p95_within_distance: float
    min_between_distance: float
    between_within_ratio: float
    pair_callable_fraction: float
    status: str


def genotype_view(genotypes: np.ndarray, model: str) -> np.ndarray:
    result = np.asarray(genotypes, dtype=np.int8)
    if model == "inbred":
        result = result.copy()
        result[result == 1] = -1
    elif model not in {"dosage", "phased"}:
        raise ValueError(f"Unknown genotype model: {model}")
    return result


def pairwise_distance_snp_per_mb(
    genotypes: np.ndarray,
    n_eligible_sites: int,
    span_bp: int,
    model: str = "inbred",
    max_sites: int = 5000,
    min_compared_sites: int = 10,
) -> DistanceResult:
    """Missing-aware genotype mismatch distance scaled to all eligible sites."""
    if genotypes.ndim != 2:
        raise ValueError("genotypes must have shape sites x samples")
    n_samples = genotypes.shape[1]
    take = stable_subsample_indices(genotypes.shape[0], max_sites)
    matrix = genotype_view(genotypes[take], model)
    if matrix.shape[0] == 0:
        distance = np.full((n_samples, n_samples), np.inf, dtype=float)
        np.fill_diagonal(distance, 0.0)
        empty = np.zeros((n_samples, n_samples), dtype=np.int64)
        return DistanceResult(distance, empty, empty, n_eligible_sites, 0, span_bp)
    valid = (matrix >= 0).astype(np.float32)
    compared_f = valid.T @ valid
    matches = np.zeros((n_samples, n_samples), dtype=np.float32)
    for state in (0, 1, 2):
        state_matrix = (matrix == state).astype(np.float32)
        matches += state_matrix.T @ state_matrix
    mismatch_f = compared_f - matches
    compared = np.rint(compared_f).astype(np.int64)
    mismatch = np.rint(mismatch_f).astype(np.int64)
    with np.errstate(divide="ignore", invalid="ignore"):
        fraction = np.divide(
            mismatch_f,
            compared_f,
            out=np.full_like(mismatch_f, np.nan, dtype=np.float32),
            where=compared_f >= min_compared_sites,
        )
    distance = fraction.astype(float) * n_eligible_sites / max(span_bp / 1_000_000.0, 1e-12)
    distance[~np.isfinite(distance)] = np.inf
    distance = np.maximum(distance, distance.T)
    np.fill_diagonal(distance, 0.0)
    return DistanceResult(distance, compared, mismatch, n_eligible_sites, matrix.shape[0], span_bp)


def _finite_values(matrix, pairs):
    values = [float(matrix[i, j]) for i, j in pairs if math.isfinite(float(matrix[i, j]))]
    return values


def complete_linkage_clusters(
    distance_result: DistanceResult,
    sample_names: list[str],
    threshold: float,
    min_hap_fraction: float,
    min_hap_samples: int,
) -> ClusterResult:
    distance = distance_result.distance
    n = len(sample_names)
    if n == 0:
        return ClusterResult([], np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32), {}, 0, 0, 0.0,
                             math.nan, math.nan, math.nan, math.nan, 0.0, "NO_SAMPLES")
    if n == 1:
        raw = np.ones(1, dtype=np.int32)
    else:
        finite = distance[np.isfinite(distance) & ~np.eye(n, dtype=bool)]
        replacement = max(threshold * 10.0 + 1.0, float(finite.max()) + threshold + 1.0 if finite.size else threshold * 10.0 + 1.0)
        work = np.where(np.isfinite(distance), distance, replacement)
        work = (work + work.T) / 2.0
        np.fill_diagonal(work, 0.0)
        condensed = squareform(work, checks=False)
        tree = linkage(condensed, method="complete", optimal_ordering=False)
        raw = fcluster(tree, t=threshold, criterion="distance").astype(np.int32)

    members = defaultdict(list)
    for index, value in enumerate(raw.tolist()):
        members[value].append(index)
    ordered = sorted(members.items(), key=lambda item: (-len(item[1]), min(sample_names[i] for i in item[1])))
    required = max(min_hap_samples, int(math.ceil(min_hap_fraction * n)))
    major_raw = [value for value, indices in ordered if len(indices) >= required]
    major_number = {value: rank + 1 for rank, value in enumerate(major_raw)}
    minor_number = {value: rank + 1 for rank, (value, indices) in enumerate((item for item in ordered if item[0] not in major_number))}
    labels = []
    numeric = np.zeros(n, dtype=np.int32)
    groups = defaultdict(list)
    for index, value in enumerate(raw.tolist()):
        if value in major_number:
            label = f"Hap{major_number[value]:02d}"
            numeric[index] = major_number[value]
        else:
            label = f"Minor{minor_number[value]:02d}"
        labels.append(label)
        groups[label].append(index)

    within_pairs = []
    between_pairs = []
    major_labels = sorted({x for x in numeric.tolist() if x > 0})
    for i in range(n):
        for j in range(i + 1, n):
            if numeric[i] > 0 and numeric[i] == numeric[j]:
                within_pairs.append((i, j))
            elif numeric[i] > 0 and numeric[j] > 0 and numeric[i] != numeric[j]:
                between_pairs.append((i, j))
    within = _finite_values(distance, within_pairs)
    between = _finite_values(distance, between_pairs)
    max_within = max(within) if within else 0.0 if major_labels else math.nan
    p95_within = float(np.quantile(within, 0.95)) if within else max_within
    min_between = min(between) if between else math.nan
    ratio = min_between / max(max_within, 1e-12) if math.isfinite(min_between) and math.isfinite(max_within) else math.nan
    assigned = int(np.sum(numeric > 0))
    tri = np.triu_indices(n, 1)
    callable_fraction = float(np.mean(np.isfinite(distance_result.distance[tri]))) if len(tri[0]) else 1.0
    diameter_ok = not within or max_within <= threshold + 1e-8
    status = "PASS" if major_labels and diameter_ok else "NO_MAJOR_HAPLOTYPE" if not major_labels else "FAIL_DIAMETER"
    return ClusterResult(
        labels, numeric, raw, dict(groups), len(major_labels), assigned, assigned / n,
        max_within, p95_within, min_between, ratio, callable_fraction, status,
    )


def association_status(labels: list[str], n_samples: int, min_frequency=0.10, min_assigned_fraction=0.50, min_haplotypes=2):
    counts = Counter(label for label in labels if label.startswith("Hap"))
    major = {label: count for label, count in counts.items() if count / max(1, n_samples) >= min_frequency}
    assigned = sum(major.values())
    reasons = []
    if len(major) < min_haplotypes:
        reasons.append("FEWER_THAN_MIN_MAJOR_HAPLOTYPES")
    if assigned / max(1, n_samples) < min_assigned_fraction:
        reasons.append("LOW_MAJOR_ASSIGNED_FRACTION")
    return {
        "association_ready": not reasons,
        "major_haplotypes": major,
        "n_major_haplotypes": len(major),
        "major_assigned_count": assigned,
        "major_assigned_fraction": assigned / max(1, n_samples),
        "reasons": reasons,
    }


def pca_for_block(genotypes: np.ndarray, positions: np.ndarray, sample_names: list[str], n_components=5, max_sites=5000):
    from sklearn.decomposition import PCA

    take = stable_subsample_indices(genotypes.shape[0], max_sites)
    if len(take) < 2 or len(sample_names) < 2:
        return [], [], []
    matrix = np.asarray(genotypes[take], dtype=float).T
    matrix[matrix < 0] = np.nan
    means = np.nanmean(matrix, axis=0)
    means[~np.isfinite(means)] = 0.0
    missing = np.where(~np.isfinite(matrix))
    matrix[missing] = means[missing[1]]
    std = np.std(matrix, axis=0)
    keep = std > 0
    if int(keep.sum()) < 2:
        return [], [], []
    matrix = (matrix[:, keep] - means[keep]) / std[keep]
    k = min(n_components, matrix.shape[0] - 1, matrix.shape[1])
    if k < 1:
        return [], [], []
    pca = PCA(n_components=k, svd_solver="randomized", random_state=1731)
    scores = pca.fit_transform(matrix)
    kept_positions = np.asarray(positions)[take][keep]
    score_rows = [
        {"SAMPLE": sample_names[i], **{f"PC{j+1}": float(scores[i, j]) for j in range(k)}}
        for i in range(len(sample_names))
    ]
    variance_rows = [{"PC": f"PC{j+1}", "EXPLAINED_VARIANCE_RATIO": float(pca.explained_variance_ratio_[j])} for j in range(k)]
    loading_rows = []
    for component in range(min(3, k)):
        order = np.argsort(np.abs(pca.components_[component]))[::-1][:min(50, pca.components_.shape[1])]
        for rank, index in enumerate(order, 1):
            loading_rows.append({
                "PC": f"PC{component+1}", "RANK": rank, "POSITION0": int(kept_positions[index]),
                "LOADING": float(pca.components_[component, index]),
            })
    return score_rows, variance_rows, loading_rows
