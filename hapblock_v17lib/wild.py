from __future__ import annotations

from collections import Counter, defaultdict

import numpy as np

from .blocks import BlockResult, interval_genotypes
from .cluster import pairwise_distance_snp_per_mb


def project_wild(
    blocks: list[BlockResult], caches, cultivar_samples, wild_samples,
    primary_maf, genotype_model, max_sites, min_sites,
):
    if not wild_samples:
        return [], []
    n_cultivar = len(cultivar_samples)
    all_indices = list(range(n_cultivar + len(wild_samples)))
    rows = []
    summaries = []
    for block in blocks:
        cache = caches[block.chrom]
        genotype, _, n_sites = interval_genotypes(
            cache, block.start0, block.end0, all_indices, primary_maf, 0.5, False
        )
        if n_sites < min_sites:
            continue
        distance = pairwise_distance_snp_per_mb(
            genotype, n_sites, block.length_bp, genotype_model, max_sites, min_sites
        ).distance
        hap_members = defaultdict(list)
        for index, label in enumerate(block.labels):
            if label.startswith("Hap"):
                hap_members[label].append(index)
        closest_counts = Counter()
        for wild_offset, sample in enumerate(wild_samples):
            wi = n_cultivar + wild_offset
            values = []
            for hap, members in sorted(hap_members.items()):
                d = [float(distance[wi, member]) for member in members if np.isfinite(distance[wi, member])]
                median = float(np.median(d)) if d else float("nan")
                mean = float(np.mean(d)) if d else float("nan")
                values.append((hap, median))
                rows.append({
                    "BLOCK_ID": block.block_id, "WILD_SAMPLE": sample, "CULTIVAR_HAPLOTYPE": hap,
                    "MEAN_DISTANCE_SNP_PER_MB": mean, "MEDIAN_DISTANCE_SNP_PER_MB": median,
                    "N_CULTIVAR_CARRIERS": len(members), "N_SITES": n_sites,
                })
            finite = [(hap, value) for hap, value in values if np.isfinite(value)]
            if finite:
                best_hap, best = min(finite, key=lambda item: (item[1], item[0]))
                closest_counts[best_hap] += 1
                for row in rows[-len(values):]:
                    row["IS_CLOSEST_HAPLOTYPE"] = "YES" if row["CULTIVAR_HAPLOTYPE"] == best_hap else "NO"
                    row["DISTANCE_DELTA_FROM_BEST"] = row["MEDIAN_DISTANCE_SNP_PER_MB"] - best if np.isfinite(row["MEDIAN_DISTANCE_SNP_PER_MB"]) else "NA"
        for hap, members in sorted(hap_members.items()):
            summaries.append({
                "BLOCK_ID": block.block_id, "CULTIVAR_HAPLOTYPE": hap,
                "N_CULTIVAR_CARRIERS": len(members), "N_WILD_CLOSEST": closest_counts[hap],
                "WILD_CLOSEST_FRACTION": closest_counts[hap] / len(wild_samples),
                "INTERPRETATION": "WILD_AFFINITY_ONLY_NOT_DIRECT_ANCESTRY",
            })
    return rows, summaries
