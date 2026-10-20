from __future__ import annotations

import math
from collections import Counter
from pathlib import Path

from .cluster import association_status
from .common import write_tsv


BLOCK_FIELDS = [
    "BLOCK_ID", "BLOCK_LEVEL", "PARENT_BLOCK_ID", "CHROM", "START0", "END", "LENGTH_BP",
    "N_ELIGIBLE_SITES", "N_SAMPLED_SITES", "SUPPORT_SCALES", "SUPPORT_MAFS",
    "N_MAJOR_HAPLOTYPES", "ASSIGNED_COUNT", "ASSIGNED_FRACTION",
    "MAX_WITHIN_DISTANCE_SNP_PER_MB", "P95_WITHIN_DISTANCE_SNP_PER_MB",
    "MIN_BETWEEN_DISTANCE_SNP_PER_MB", "BETWEEN_WITHIN_RATIO", "PAIR_CALLABLE_FRACTION",
    "LOCAL_PASS_FRACTION", "N_INTERNAL_GAPS", "CLUSTER_STATUS",
    "ASSOCIATION_READY", "ASSOCIATION_N_MAJOR_HAPS", "ASSOCIATION_ASSIGNED_FRACTION",
    "ASSOCIATION_FAILURE_REASONS",
]


def _number(value):
    return value if isinstance(value, (int, str)) or math.isfinite(value) else "NA"


def block_rows(blocks, n_samples, association_frequency, association_assigned, association_min_haps):
    rows = []
    association = {}
    for block in blocks:
        assoc = association_status(block.labels, n_samples, association_frequency, association_assigned, association_min_haps)
        association[block.block_id] = assoc
        rows.append({
            "BLOCK_ID": block.block_id, "BLOCK_LEVEL": block.level, "PARENT_BLOCK_ID": block.parent_id,
            "CHROM": block.chrom, "START0": block.start0, "END": block.end0, "LENGTH_BP": block.length_bp,
            "N_ELIGIBLE_SITES": block.n_sites, "N_SAMPLED_SITES": block.n_sampled_sites,
            "SUPPORT_SCALES": ",".join(map(str, block.support_scales)),
            "SUPPORT_MAFS": ",".join(f"{x:g}" for x in block.support_mafs),
            "N_MAJOR_HAPLOTYPES": block.cluster.n_major_haps,
            "ASSIGNED_COUNT": block.cluster.assigned_count, "ASSIGNED_FRACTION": block.cluster.assigned_fraction,
            "MAX_WITHIN_DISTANCE_SNP_PER_MB": _number(block.cluster.max_within_distance),
            "P95_WITHIN_DISTANCE_SNP_PER_MB": _number(block.cluster.p95_within_distance),
            "MIN_BETWEEN_DISTANCE_SNP_PER_MB": _number(block.cluster.min_between_distance),
            "BETWEEN_WITHIN_RATIO": _number(block.cluster.between_within_ratio),
            "PAIR_CALLABLE_FRACTION": block.cluster.pair_callable_fraction,
            "LOCAL_PASS_FRACTION": block.local_pass_fraction, "N_INTERNAL_GAPS": len(block.gaps),
            "CLUSTER_STATUS": block.cluster.status,
            "ASSOCIATION_READY": "YES" if assoc["association_ready"] else "NO",
            "ASSOCIATION_N_MAJOR_HAPS": assoc["n_major_haplotypes"],
            "ASSOCIATION_ASSIGNED_FRACTION": assoc["major_assigned_fraction"],
            "ASSOCIATION_FAILURE_REASONS": ",".join(assoc["reasons"]) if assoc["reasons"] else "NONE",
        })
    return rows, association


def write_core_outputs(outdir, prefix, blocks, samples, association_frequency, association_assigned, association_min_haps):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    rows, association = block_rows(blocks, len(samples), association_frequency, association_assigned, association_min_haps)
    write_tsv(out / f"{prefix}.blocks.tsv", BLOCK_FIELDS, rows)
    with (out / f"{prefix}.blocks.bed").open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(f"{row['CHROM']}\t{row['START0']}\t{row['END']}\t{row['BLOCK_ID']}\t{row['BLOCK_LEVEL']}\t{row['ASSOCIATION_READY']}\n")
    hierarchy = [{
        "BLOCK_ID": block.block_id, "PARENT_BLOCK_ID": block.parent_id, "BLOCK_LEVEL": block.level,
        "CHROM": block.chrom, "START0": block.start0, "END": block.end0,
    } for block in blocks]
    write_tsv(out / f"{prefix}.block_hierarchy.tsv", ["BLOCK_ID", "PARENT_BLOCK_ID", "BLOCK_LEVEL", "CHROM", "START0", "END"], hierarchy)
    gap_rows = []
    for block in blocks:
        for start0, end0, reason in block.gaps:
            gap_rows.append({"BLOCK_ID": block.block_id, "CHROM": block.chrom, "START0": start0, "END": end0, "LENGTH_BP": end0-start0, "REASON": reason})
    write_tsv(out / f"{prefix}.internal_gaps.tsv", ["BLOCK_ID", "CHROM", "START0", "END", "LENGTH_BP", "REASON"], gap_rows)

    freq_rows, sample_rows = [], []
    for block in blocks:
        counts = Counter(block.labels)
        assoc_major = association[block.block_id]["major_haplotypes"]
        for label, count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
            freq_rows.append({
                "BLOCK_ID": block.block_id, "HAPLOTYPE": label, "N_SAMPLES": count,
                "FREQUENCY_ALL_CULTIVARS": count / len(samples),
                "IS_ASSOCIATION_MAJOR_HAPLOTYPE": "YES" if label in assoc_major else "NO",
            })
        for sample, label in zip(samples, block.labels):
            sample_rows.append({"SAMPLE": sample, "BLOCK_ID": block.block_id, "HAPLOTYPE": label})
    write_tsv(out / f"{prefix}.block_haplotype_frequency.tsv",
              ["BLOCK_ID", "HAPLOTYPE", "N_SAMPLES", "FREQUENCY_ALL_CULTIVARS", "IS_ASSOCIATION_MAJOR_HAPLOTYPE"], freq_rows)
    write_tsv(out / f"{prefix}.block_haplotype_samples.tsv", ["SAMPLE", "BLOCK_ID", "HAPLOTYPE"], sample_rows)

    with (out / f"{prefix}.sample_by_block_haplotype.tsv").open("w", encoding="utf-8") as handle:
        handle.write("SAMPLE\t" + "\t".join(block.block_id for block in blocks) + "\n")
        for i, sample in enumerate(samples):
            handle.write(sample + "\t" + "\t".join(block.labels[i] for block in blocks) + "\n")

    ready = [block for block in blocks if association[block.block_id]["association_ready"]]
    ready_ids = {block.block_id for block in ready}
    write_tsv(out / f"{prefix}.association_ready_blocks.tsv", BLOCK_FIELDS, [row for row in rows if row["BLOCK_ID"] in ready_ids])
    write_tsv(out / f"{prefix}.association_excluded_blocks.tsv", BLOCK_FIELDS, [row for row in rows if row["BLOCK_ID"] not in ready_ids])
    with (out / f"{prefix}.association_block_factor_matrix.tsv").open("w", encoding="utf-8") as handle:
        handle.write("SAMPLE\t" + "\t".join(block.block_id for block in ready) + "\n")
        for i, sample in enumerate(samples):
            handle.write(sample + "\t" + "\t".join(block.labels[i] for block in ready) + "\n")
    onehot_columns = []
    for block in ready:
        for hap in sorted(association[block.block_id]["major_haplotypes"]):
            onehot_columns.append((block, hap, f"{block.block_id}__{hap}"))
    with (out / f"{prefix}.association_haplotype_onehot_matrix.tsv").open("w", encoding="utf-8") as handle:
        handle.write("SAMPLE\t" + "\t".join(column for _, _, column in onehot_columns) + "\n")
        for i, sample in enumerate(samples):
            handle.write(sample + "\t" + "\t".join("1" if block.labels[i] == hap else "0" for block, hap, _ in onehot_columns) + "\n")
    return rows, association


def write_report(path, prefix, mode, manifest, block_rows, association, wild_enabled):
    coarse = [row for row in block_rows if row["BLOCK_LEVEL"] == "COARSE"]
    fine = [row for row in block_rows if row["BLOCK_LEVEL"] == "FINE"]
    ready = [row for row in block_rows if row["ASSOCIATION_READY"] == "YES"]
    text = [
        f"# HapBlock V17 report: {prefix}", "",
        f"- Mode: `{mode}`",
        f"- Cultivar samples: {manifest['n_cultivar']}",
        f"- Wild samples: {manifest['n_wild']}",
        f"- Retained biallelic SNP records: {manifest['n_retained_biallelic_snps']}",
        f"- Coarse blocks: {len(coarse)}",
        f"- Fine blocks: {len(fine)}",
        f"- Association-ready blocks: {len(ready)}", "",
        "A block is association-ready only when it contains at least the configured number of major haplotypes, each major haplotype reaches the configured frequency, and their combined assigned fraction reaches the configured threshold.", "",
        "PCA is a QC and visualization layer. Final haplotype assignments are defined by complete-linkage clustering and verified against the original SNP-distance threshold.", "",
    ]
    if not wild_enabled:
        text += ["Wild samples were not provided. No ancestral/derived interpretation was attempted.", ""]
    text += ["## Association-ready blocks", ""]
    if ready:
        for row in ready:
            text.append(f"- `{row['BLOCK_ID']}`: {row['CHROM']}:{int(row['START0'])+1:,}-{int(row['END']):,}; {row['ASSOCIATION_N_MAJOR_HAPS']} major haplotypes; assigned fraction {float(row['ASSOCIATION_ASSIGNED_FRACTION']):.3f}")
    else:
        text.append("- None")
    Path(path).write_text("\n".join(text) + "\n", encoding="utf-8")


def write_recent_outputs(outdir, prefix, samples, block_ids, assignments, rows, frequencies,
                         min_frequency, min_assigned_fraction, min_haplotypes):
    out = Path(outdir)
    write_tsv(out / f"{prefix}.backbone_to_recent_subhap.tsv",
              ["SAMPLE", "BLOCK_ID", "BACKBONE_HAPLOTYPE", "RECENT_SUBHAPLOTYPE", "SUBHAP_STATUS"], rows)
    write_tsv(out / f"{prefix}.recent_subhaplotype_frequency.tsv",
              ["BLOCK_ID", "BACKBONE_HAPLOTYPE", "RECENT_SUBHAPLOTYPE", "N_SAMPLES",
               "FREQUENCY_WITHIN_BACKBONE", "FREQUENCY_ALL_CULTIVARS", "N_POLYMORPHIC_SITES", "STATUS"], frequencies)
    with (out / f"{prefix}.sample_by_recent_subhaplotype.tsv").open("w", encoding="utf-8") as handle:
        handle.write("SAMPLE\t" + "\t".join(block_ids) + "\n")
        for index, sample in enumerate(samples):
            handle.write(sample + "\t" + "\t".join(assignments.get(block_id, ["NA"] * len(samples))[index] for block_id in block_ids) + "\n")
    ready_rows = []
    columns = []
    for block_id in block_ids:
        labels = assignments.get(block_id, [])
        counts = Counter(
            label for label in labels
            if label not in {"NA", ""} and ".MinorSubHap" not in label and ".UnresolvedSubHap" not in label
        )
        major = {label: count for label, count in counts.items() if count / max(1, len(samples)) >= min_frequency}
        assigned = sum(major.values())
        ready = len(major) >= min_haplotypes and assigned / max(1, len(samples)) >= min_assigned_fraction
        ready_rows.append({
            "BLOCK_ID": block_id, "N_MAJOR_RECENT_SUBHAPS": len(major),
            "MAJOR_ASSIGNED_COUNT": assigned, "MAJOR_ASSIGNED_FRACTION": assigned / max(1, len(samples)),
            "ASSOCIATION_READY": "YES" if ready else "NO",
        })
        if ready:
            for label in sorted(major):
                columns.append((block_id, label, f"{block_id}__{label}"))
    write_tsv(out / f"{prefix}.recent_subhap_association_status.tsv",
              ["BLOCK_ID", "N_MAJOR_RECENT_SUBHAPS", "MAJOR_ASSIGNED_COUNT", "MAJOR_ASSIGNED_FRACTION", "ASSOCIATION_READY"], ready_rows)
    with (out / f"{prefix}.association_recent_subhap_onehot_matrix.tsv").open("w", encoding="utf-8") as handle:
        handle.write("SAMPLE\t" + "\t".join(column for _, _, column in columns) + "\n")
        for index, sample in enumerate(samples):
            values = []
            for block_id, label, _ in columns:
                values.append("1" if assignments[block_id][index] == label else "0")
            handle.write(sample + "\t" + "\t".join(values) + "\n")
