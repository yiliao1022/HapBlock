from __future__ import annotations

from collections import Counter

import numpy as np

from .blocks import BlockResult, call_interval


def call_recent_subhaps(
    blocks: list[BlockResult],
    caches,
    cultivar_samples,
    threshold,
    genotype_model,
    min_sites,
    max_sites,
    min_subhap_fraction,
    min_subhap_samples,
):
    assignments = {}
    rows = []
    frequencies = []
    for block in blocks:
        cache = caches[block.chrom]
        combined = ["NA"] * len(cultivar_samples)
        for backbone in sorted({x for x in block.labels if x.startswith("Hap")}):
            member_indices = [i for i, value in enumerate(block.labels) if value == backbone]
            if not member_indices:
                continue
            member_names = [cultivar_samples[i] for i in member_indices]
            cluster, genotype, _ = call_interval(
                cache, block.start0, block.end0, member_indices, member_names,
                1e-12, threshold, genotype_model, min_sites, max_sites,
                min_subhap_fraction, min_subhap_samples, maf_max=0.5, recalc_subset_maf=True,
            )
            if cluster is None:
                local = ["UnresolvedSubHap"] * len(member_indices)
                status = "LOW_INFORMATION_UNRESOLVED"
                n_sites = int(genotype.shape[0])
            else:
                local = [value.replace("Hap", "SubHap", 1) if value.startswith("Hap") else value.replace("Minor", "MinorSubHap", 1) for value in cluster.labels]
                status = cluster.status
                n_sites = int(genotype.shape[0])
            for local_index, global_index in enumerate(member_indices):
                combined[global_index] = f"{backbone}.{local[local_index]}"
                rows.append({
                    "SAMPLE": cultivar_samples[global_index], "BLOCK_ID": block.block_id,
                    "BACKBONE_HAPLOTYPE": backbone, "RECENT_SUBHAPLOTYPE": combined[global_index],
                    "SUBHAP_STATUS": status,
                })
            counts = Counter(local)
            for label, count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
                frequencies.append({
                    "BLOCK_ID": block.block_id, "BACKBONE_HAPLOTYPE": backbone,
                    "RECENT_SUBHAPLOTYPE": f"{backbone}.{label}", "N_SAMPLES": count,
                    "FREQUENCY_WITHIN_BACKBONE": count / len(member_indices),
                    "FREQUENCY_ALL_CULTIVARS": count / len(cultivar_samples),
                    "N_POLYMORPHIC_SITES": n_sites, "STATUS": status,
                })
        assignments[block.block_id] = combined
    return assignments, rows, frequencies
