#!/usr/bin/env python3
"""HapBlock V17: hierarchical cultivar haplotype blocks and association-ready groups."""
from __future__ import annotations

import argparse
import json
import math
import os
import shlex
import shutil
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from hapblock_v17lib import __version__
from hapblock_v17lib.blocks import (
    fine_children_for_parent, interval_genotypes, merge_scale_blocks, merge_window_calls,
    scan_windows, validate_candidates,
)
from hapblock_v17lib.cache import build_cache, open_cache
from hapblock_v17lib.cluster import genotype_view, pca_for_block
from hapblock_v17lib.common import chrom_key, parse_float_csv, parse_int_csv, write_json, write_tsv
from hapblock_v17lib.outputs import write_core_outputs, write_recent_outputs, write_report
from hapblock_v17lib.recent import call_recent_subhaps
from hapblock_v17lib.wild import project_wild


def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--version", action="version", version=f"HapBlock V{__version__}")
    p.add_argument("--vcf", required=True, help="Unfiltered, coordinate-sorted diploid SNP VCF/VCF.gz")
    p.add_argument("--cultivar-list", required=True)
    p.add_argument("--wild-list", default="")
    p.add_argument("--outdir", required=True)
    p.add_argument("--prefix", default="hapblock_v17")
    p.add_argument("--genotype-model", choices=("inbred", "dosage"), default="inbred")
    p.add_argument("--backbone-mafs", default="0.05,0.10")
    p.add_argument("--coarse-windows", default="1000000,500000")
    p.add_argument("--fine-windows", default="100000,50000,25000,10000")
    p.add_argument("--fine-context-bp", type=int, default=100000)
    p.add_argument("--fine-max-context-bp", type=int, default=200000)
    p.add_argument("--backbone-threshold-snp-per-mb", type=float, default=100.0)
    p.add_argument("--recent-threshold-snp-per-mb", type=float, default=20.0)
    p.add_argument("--min-window-snps", type=int, default=10)
    p.add_argument("--max-window-snps", type=int, default=2000)
    p.add_argument("--max-block-snps", type=int, default=10000)
    p.add_argument("--min-structural-hap-frequency", type=float, default=0.05)
    p.add_argument("--min-structural-hap-samples", type=int, default=3)
    p.add_argument("--min-recent-subhap-frequency", type=float, default=0.05)
    p.add_argument("--min-recent-subhap-samples", type=int, default=3)
    p.add_argument("--membership-ari", type=float, default=0.85)
    p.add_argument("--min-scale-support", type=int, default=2)
    p.add_argument("--min-coarse-block-bp", type=int, default=500000)
    p.add_argument("--min-fine-block-bp", type=int, default=10000)
    p.add_argument("--coarse-max-gap-bp", type=int, default=1000000)
    p.add_argument("--fine-max-gap-bp", type=int, default=50000)
    p.add_argument("--max-gap-fraction", type=float, default=0.05)
    p.add_argument("--min-association-hap-frequency", type=float, default=0.10)
    p.add_argument("--min-association-assigned-fraction", type=float, default=0.50)
    p.add_argument("--min-association-haplotypes", type=int, default=2)
    p.add_argument("--pca-components", type=int, default=5)
    p.add_argument("--pca-max-blocks", type=int, default=100)
    p.add_argument("--pca-max-snps", type=int, default=5000)
    p.add_argument("--skip-pca", action="store_true")
    p.add_argument("--skip-recent-subhaps", action="store_true")
    p.add_argument("--skip-wild-projection", action="store_true")
    p.add_argument("--keep-chr", default="", help="Comma-separated chromosomes to analyze after caching")
    p.add_argument("--cache-dir", default="", help="Optional persistent genotype cache directory")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--force", action="store_true", help="Replace an existing output directory/cache")
    return p.parse_args()


def validate_args(args):
    for path, label in ((args.vcf, "VCF"), (args.cultivar_list, "cultivar list")):
        if not Path(path).exists():
            raise SystemExit(f"{label} not found: {path}")
    if args.wild_list and not Path(args.wild_list).exists():
        raise SystemExit(f"wild list not found: {args.wild_list}")
    for value, label in (
        (args.min_structural_hap_frequency, "min structural hap frequency"),
        (args.min_recent_subhap_frequency, "min recent subhap frequency"),
        (args.min_association_hap_frequency, "min association hap frequency"),
        (args.min_association_assigned_fraction, "min association assigned fraction"),
        (args.max_gap_fraction, "max gap fraction"),
    ):
        if not 0 <= value <= 1:
            raise SystemExit(f"{label} must be between 0 and 1")
    if args.fine_max_context_bp < args.fine_context_bp:
        raise SystemExit("fine-max-context-bp must be >= fine-context-bp")


def window_row(call):
    cluster = call.cluster
    return {
        "CHROM": call.chrom, "START0": call.start0, "END": call.end0,
        "CONTEXT_START0": call.context_start0, "CONTEXT_END": call.context_end0,
        "WINDOW_BP": call.window_bp, "STEP_BP": call.step_bp, "MAF_MIN": call.maf_min,
        "N_SITES": call.n_sites, "N_SAMPLED_SITES": call.n_sampled_sites, "STATUS": call.status,
        "N_MAJOR_HAPLOTYPES": cluster.n_major_haps if cluster else 0,
        "ASSIGNED_FRACTION": cluster.assigned_fraction if cluster else 0,
        "MAX_WITHIN_DISTANCE_SNP_PER_MB": cluster.max_within_distance if cluster and math.isfinite(cluster.max_within_distance) else "NA",
        "MIN_BETWEEN_DISTANCE_SNP_PER_MB": cluster.min_between_distance if cluster and math.isfinite(cluster.min_between_distance) else "NA",
    }


def pca_outputs(out, prefix, blocks, caches, samples, primary_maf, model, n_components, max_blocks, max_sites):
    score_rows, variance_rows, loading_rows = [], [], []
    sample_index = {sample: index for index, sample in enumerate(samples)}
    selected = sorted(blocks, key=lambda block: (-block.length_bp, block.block_id))[:max_blocks]
    pdf = None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        pdf = PdfPages(out / f"{prefix}.block_pca_plots.pdf")
    except Exception:
        plt = None
    for block in selected:
        genotype, positions, _ = interval_genotypes(
            caches[block.chrom], block.start0, block.end0, list(range(len(samples))), primary_maf, 0.5, False
        )
        genotype = genotype_view(genotype, model)
        scores, variance, loadings = pca_for_block(genotype, positions, samples, n_components, max_sites)
        for row in scores:
            row.update({"BLOCK_ID": block.block_id, "HAPLOTYPE": block.labels[sample_index[row["SAMPLE"]]]})
        for row in variance:
            row["BLOCK_ID"] = block.block_id
        for row in loadings:
            row.update({"BLOCK_ID": block.block_id, "CHROM": block.chrom})
        score_rows.extend(scores); variance_rows.extend(variance); loading_rows.extend(loadings)
        if pdf is not None and scores and "PC2" in scores[0]:
            fig, ax = plt.subplots(figsize=(7, 6))
            labels = sorted(set(row["HAPLOTYPE"] for row in scores))
            cmap = plt.get_cmap("tab20")
            for i, label in enumerate(labels):
                subset = [row for row in scores if row["HAPLOTYPE"] == label]
                ax.scatter([row["PC1"] for row in subset], [row["PC2"] for row in subset], s=14, alpha=.75, label=label, color=cmap(i % 20))
            ax.set_title(f"{block.block_id} | {block.chrom}:{block.start0+1}-{block.end0}")
            ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.legend(fontsize=6, ncol=2)
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    if pdf is not None:
        pdf.close()
    score_fields = ["BLOCK_ID", "SAMPLE", "HAPLOTYPE"] + [f"PC{i+1}" for i in range(n_components)]
    write_tsv(out / f"{prefix}.block_pca_scores.tsv", score_fields, score_rows)
    write_tsv(out / f"{prefix}.block_pca_variance.tsv", ["BLOCK_ID", "PC", "EXPLAINED_VARIANCE_RATIO"], variance_rows)
    write_tsv(out / f"{prefix}.block_pca_top_loadings.tsv", ["BLOCK_ID", "CHROM", "PC", "RANK", "POSITION0", "LOADING"], loading_rows)


def plot_landscape(out, prefix, blocks, association):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    chromosomes = sorted({block.chrom for block in blocks}, key=chrom_key)
    if not chromosomes:
        return
    fig, axes = plt.subplots(len(chromosomes), 1, figsize=(14, max(3, len(chromosomes) * 1.5)), squeeze=False)
    for ax, chrom in zip(axes[:, 0], chromosomes):
        subset = [block for block in blocks if block.chrom == chrom]
        for block in subset:
            y = 1 if block.level == "COARSE" else 0
            color = "#3a86ff" if association[block.block_id]["association_ready"] else "#9aa0a6"
            ax.plot([block.start0 / 1e6, block.end0 / 1e6], [y, y], lw=7 if y else 4, color=color, solid_capstyle="butt")
        ax.set_yticks([0, 1], ["Fine", "Coarse"]); ax.set_ylabel(chrom); ax.grid(axis="x", alpha=.2)
    axes[-1, 0].set_xlabel("Position (Mb)")
    fig.suptitle(f"{prefix} hierarchical haplotype blocks (blue = association-ready)")
    fig.tight_layout(); fig.savefig(out / f"{prefix}.block_landscape.pdf"); fig.savefig(out / f"{prefix}.block_landscape.png", dpi=180); plt.close(fig)


def main():
    args = parse_args(); validate_args(args)
    out = Path(args.outdir).resolve()
    if out.exists() and any(out.iterdir()) and not args.resume:
        if not args.force:
            raise SystemExit(f"Output directory is not empty: {out}. Use --resume to reuse its validated cache or --force to replace it")
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    run_info = out / "00_run_info"; run_info.mkdir(exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(run_info / ".matplotlib")
    os.environ.setdefault("XDG_CACHE_HOME", str(run_info / ".cache"))
    write_json(run_info / "parameters.json", {**vars(args), "version": __version__})
    (run_info / "command.txt").write_text(shlex.join([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]]) + "\n", encoding="utf-8")
    cache_dir = Path(args.cache_dir).resolve() if args.cache_dir else out / "01_genotype_cache"
    started = time.time()
    print("[01] building/reusing genotype cache", file=sys.stderr, flush=True)
    build_cache(args.vcf, args.cultivar_list, args.wild_list, cache_dir, args.resume, args.force)
    manifest, caches = open_cache(cache_dir)
    write_tsv(
        run_info / "sample_qc.tsv",
        ["SAMPLE", "GROUP", "CALLED_GENOTYPES", "MISSING_GENOTYPES", "CALL_RATE", "HETEROZYGOUS_GENOTYPES", "HETEROZYGOSITY_AMONG_CALLED"],
        [
            {
                "SAMPLE": row["sample"], "GROUP": row["group"],
                "CALLED_GENOTYPES": row["called_genotypes"], "MISSING_GENOTYPES": row["missing_genotypes"],
                "CALL_RATE": row["call_rate"], "HETEROZYGOUS_GENOTYPES": row["heterozygous_genotypes"],
                "HETEROZYGOSITY_AMONG_CALLED": row["heterozygosity_among_called"],
            }
            for row in manifest.get("sample_qc", [])
        ],
    )
    samples = manifest["cultivar_samples"]
    mode = "CULTIVAR_WILD" if manifest["n_wild"] else "CULTIVAR_ONLY"
    chromosomes = sorted(caches, key=chrom_key)
    if args.keep_chr:
        wanted_chroms = {x.strip() for x in args.keep_chr.split(",") if x.strip()}
        absent = wanted_chroms - set(chromosomes)
        if absent:
            raise SystemExit("Requested chromosomes absent from cache: " + ",".join(sorted(absent)))
        chromosomes = [chrom for chrom in chromosomes if chrom in wanted_chroms]
    backbone_mafs = sorted(set(parse_float_csv(args.backbone_mafs)))
    if any(x < 0 or x > .5 for x in backbone_mafs):
        raise SystemExit("backbone MAF values must be between 0 and 0.5")
    primary_maf = backbone_mafs[0]
    coarse_windows = sorted(set(parse_int_csv(args.coarse_windows)), reverse=True)
    fine_windows = sorted(set(parse_int_csv(args.fine_windows)), reverse=True)
    cultivar_indices = list(range(len(samples)))
    all_window_rows = []
    coarse_scale_blocks = []
    print(f"[02-04] mode={mode} chromosomes={len(chromosomes)}", file=sys.stderr, flush=True)
    for chrom in chromosomes:
        cache = caches[chrom]
        for maf in backbone_mafs:
            for window in coarse_windows:
                step = max(1, window // 2)
                calls = scan_windows(
                    cache, cultivar_indices, samples, window, step, maf,
                    args.backbone_threshold_snp_per_mb, args.genotype_model,
                    args.min_window_snps, args.max_window_snps,
                    args.min_structural_hap_frequency, args.min_structural_hap_samples,
                    args.fine_context_bp, args.fine_max_context_bp,
                )
                all_window_rows.extend(window_row(call) for call in calls)
                coarse_scale_blocks.extend(merge_window_calls(
                    calls, args.membership_ari, args.coarse_max_gap_bp,
                    args.max_gap_fraction, args.min_coarse_block_bp,
                ))
                print(f"[coarse] {chrom} maf={maf:g} window={window} calls={len(calls)}", file=sys.stderr, flush=True)
    candidates = merge_scale_blocks(
        coarse_scale_blocks, args.membership_ari, args.coarse_max_gap_bp, args.max_gap_fraction
    )
    parents = validate_candidates(
        candidates, caches, cultivar_indices, samples, "COARSE", "CB", primary_maf,
        args.backbone_threshold_snp_per_mb, args.genotype_model, args.min_window_snps,
        args.max_block_snps, args.min_structural_hap_frequency, args.min_structural_hap_samples,
        args.min_scale_support,
    )
    print(f"[04] validated coarse blocks={len(parents)}", file=sys.stderr, flush=True)

    children = []
    print("[05] fine blocks and 10-kb boundary refinement", file=sys.stderr, flush=True)
    for parent in parents:
        fine_by_scale = {}
        cache = caches[parent.chrom]
        for window in fine_windows:
            step = 10000 if window <= 25000 else max(1, window // 2)
            calls = scan_windows(
                cache, cultivar_indices, samples, window, step, primary_maf,
                args.backbone_threshold_snp_per_mb, args.genotype_model,
                args.min_window_snps, args.max_window_snps,
                args.min_structural_hap_frequency, args.min_structural_hap_samples,
                args.fine_context_bp, args.fine_max_context_bp, (parent.start0, parent.end0),
            )
            all_window_rows.extend(window_row(call) for call in calls)
            fine_by_scale[window] = merge_window_calls(
                calls, args.membership_ari, args.fine_max_gap_bp,
                args.max_gap_fraction, args.min_fine_block_bp,
            )
        children.extend(fine_children_for_parent(
            parent, fine_by_scale, cache, cultivar_indices, samples, primary_maf,
            args.backbone_threshold_snp_per_mb, args.genotype_model, args.min_window_snps,
            args.max_block_snps, args.min_structural_hap_frequency, args.min_structural_hap_samples,
        ))
        print(f"[fine] {parent.block_id} children={sum(child.parent_id == parent.block_id for child in children)}", file=sys.stderr, flush=True)
    blocks = parents + children
    blocks.sort(key=lambda block: (chrom_key(block.chrom), block.start0, 0 if block.level == "COARSE" else 1, block.end0))
    window_fields = [
        "CHROM", "START0", "END", "CONTEXT_START0", "CONTEXT_END", "WINDOW_BP", "STEP_BP", "MAF_MIN",
        "N_SITES", "N_SAMPLED_SITES", "STATUS", "N_MAJOR_HAPLOTYPES", "ASSIGNED_FRACTION",
        "MAX_WITHIN_DISTANCE_SNP_PER_MB", "MIN_BETWEEN_DISTANCE_SNP_PER_MB",
    ]
    write_tsv(out / f"{args.prefix}.window_calls.tsv", window_fields, all_window_rows)
    block_rows, association = write_core_outputs(
        out, args.prefix, blocks, samples, args.min_association_hap_frequency,
        args.min_association_assigned_fraction, args.min_association_haplotypes,
    )

    if not args.skip_recent_subhaps and blocks:
        print("[07] recent subhaplotypes", file=sys.stderr, flush=True)
        recent_assignments, recent_rows, recent_frequencies = call_recent_subhaps(
            blocks, caches, samples, args.recent_threshold_snp_per_mb, args.genotype_model,
            args.min_window_snps, args.max_block_snps, args.min_recent_subhap_frequency,
            args.min_recent_subhap_samples,
        )
        write_recent_outputs(
            out, args.prefix, samples, [block.block_id for block in blocks], recent_assignments,
            recent_rows, recent_frequencies, args.min_association_hap_frequency,
            args.min_association_assigned_fraction, args.min_association_haplotypes,
        )

    if not args.skip_pca and blocks:
        print("[08] PCA QC", file=sys.stderr, flush=True)
        pca_outputs(out, args.prefix, blocks, caches, samples, primary_maf, args.genotype_model,
                    args.pca_components, args.pca_max_blocks, args.pca_max_snps)

    if manifest["n_wild"] and not args.skip_wild_projection and blocks:
        print("[09] projecting wild samples", file=sys.stderr, flush=True)
        projection, wild_summary = project_wild(
            blocks, caches, samples, manifest["wild_samples"], primary_maf,
            args.genotype_model, args.max_block_snps, args.min_window_snps,
        )
        write_tsv(out / f"{args.prefix}.wild_to_cultivar_haplotype.tsv",
                  ["BLOCK_ID", "WILD_SAMPLE", "CULTIVAR_HAPLOTYPE", "MEAN_DISTANCE_SNP_PER_MB",
                   "MEDIAN_DISTANCE_SNP_PER_MB", "N_CULTIVAR_CARRIERS", "N_SITES",
                   "IS_CLOSEST_HAPLOTYPE", "DISTANCE_DELTA_FROM_BEST"], projection)
        write_tsv(out / f"{args.prefix}.wild_support_summary.tsv",
                  ["BLOCK_ID", "CULTIVAR_HAPLOTYPE", "N_CULTIVAR_CARRIERS", "N_WILD_CLOSEST",
                   "WILD_CLOSEST_FRACTION", "INTERPRETATION"], wild_summary)

    plot_landscape(out, args.prefix, blocks, association)
    write_report(out / f"{args.prefix}.report.md", args.prefix, mode, manifest, block_rows, association, bool(manifest["n_wild"]))
    expected = [
        ("BLOCKS", out / f"{args.prefix}.blocks.tsv"),
        ("HAPLOTYPE_SAMPLES", out / f"{args.prefix}.block_haplotype_samples.tsv"),
        ("ASSOCIATION_BLOCKS", out / f"{args.prefix}.association_ready_blocks.tsv"),
        ("ASSOCIATION_FACTORS", out / f"{args.prefix}.association_block_factor_matrix.tsv"),
        ("ASSOCIATION_ONEHOT", out / f"{args.prefix}.association_haplotype_onehot_matrix.tsv"),
        ("RECENT_SUBHAPS", out / f"{args.prefix}.backbone_to_recent_subhap.tsv"),
        ("PCA_SCORES", out / f"{args.prefix}.block_pca_scores.tsv"),
        ("WILD_PROJECTION", out / f"{args.prefix}.wild_to_cultivar_haplotype.tsv"),
        ("REPORT", out / f"{args.prefix}.report.md"),
    ]
    write_tsv(out / f"{args.prefix}.output_index.tsv", ["OUTPUT", "PATH", "STATUS"], [
        {"OUTPUT": label, "PATH": str(path), "STATUS": "PRESENT" if path.exists() else "NOT_GENERATED"}
        for label, path in expected
    ])
    print(f"[V17 complete] blocks={len(blocks)} elapsed={time.time()-started:.1f}s output={out}", file=sys.stderr)


if __name__ == "__main__":
    main()
