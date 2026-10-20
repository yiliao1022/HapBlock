# HapBlock V17 User Manual

Version: 17.0.0  
Date: 2026-10-20

## 1. Purpose

HapBlock V17 detects cultivar haplotype blocks from a population VCF and reports:

1. the genomic position and length of each block;
2. the number and frequency of major haplotype groups in each block;
3. the samples carrying every haplotype;
4. fine child blocks within multi-megabase coarse parent blocks;
5. blocks that have enough samples for phenotype association analysis;
6. when wild samples are supplied, the cultivar haplotype closest to each wild sample.

The pipeline automatically selects one of two modes:

- `CULTIVAR_ONLY`: a cultivar list is supplied; the pipeline calls blocks, haplotypes, recent subhaplotypes, and association matrices.
- `CULTIVAR_WILD`: cultivar and wild lists are supplied. Cultivars still define every block and cultivar group; wild samples are projected only after those calls are fixed.

Names such as `Hap01` denote groups with similar diploid genotype patterns across an interval. In highly inbred material, these groups are often a useful approximation to haplotypes. V17 does not split an unphased heterozygous genotype into two chromosome-level haplotypes, so its output must not be interpreted as experimentally phased chromosome sequences.

## 2. What changed in V17

V17 introduces the following changes relative to the previous workflow:

- complete-linkage clustering controls the farthest pair of samples within a group;
- A–B–C chaining cannot join A and C when their distance exceeds the threshold;
- there is no fixed maximum number of reported haplotypes;
- the VCF is read once into a compact chromosome-level cache;
- coarse parent blocks and fine child blocks are both retained, with scans down to 10 kb;
- short low-information interruptions are recorded as internal gaps instead of automatically breaking a long block;
- common SNPs define older backbone haplotypes, followed by recent subhaplotypes from variants polymorphic within each backbone group;
- factor and one-hot matrices are written directly for phenotype association analysis.

## 3. Installation

Conda is recommended:

```bash
tar -xzf hapblock_v17.0.0.tar.gz
cd hapblock_v17_release
conda env create -f environment.yml
conda activate hapblock-v17
python3 hapblock_v17.py --version
```

Alternatively, use a Python virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 hapblock_v17.py --version
```

Python 3.10 or later is required. The core dependencies are NumPy, SciPy, scikit-learn, and Matplotlib.

## 4. Input files

### 4.1 VCF

The input must be a coordinate-sorted VCF or `.vcf.gz` file. Recommended properties are:

- an unfiltered population VCF, with no prior MAF filter;
- a `GT` FORMAT field;
- cultivar and optional wild samples in the same VCF;
- `##contig=<ID=...,length=...>` header records;
- diploid genotypes for the current `inbred` and `dosage` models.

V17 retains biallelic SNPs. It skips indels, multiallelic records, and records without `GT`. Phasing is not required.

If contig lengths are absent from the header, the position of the final retained SNP becomes the scanning endpoint. Keeping contig length records is therefore recommended.

### 4.2 Sample lists

`cultivar.list` and `wild.list` contain one VCF sample ID per line:

```text
sample_001
sample_002
sample_003
```

Blank lines and lines beginning with `#` are ignored. Duplicate IDs are removed. The lists must not overlap, and every requested ID must occur in the VCF header.

`--wild-list` is optional. If it is omitted, the run automatically uses `CULTIVAR_ONLY` mode.

## 5. Basic commands

### 5.1 Cultivars only

```bash
python3 hapblock_v17.py \
  --vcf population.unfiltered.vcf.gz \
  --cultivar-list cultivar.list \
  --outdir results_v17 \
  --prefix CA59 \
  2>&1 | tee CA59.v17.log
```

### 5.2 Cultivars plus wild samples

```bash
python3 hapblock_v17.py \
  --vcf population.unfiltered.vcf.gz \
  --cultivar-list cultivar.list \
  --wild-list wild.list \
  --outdir results_v17 \
  --prefix CA59 \
  2>&1 | tee CA59.v17.log
```

Wild samples do not contribute to cultivar block boundaries or cultivar clustering.

### 5.3 Test one chromosome first

```bash
python3 hapblock_v17.py \
  --vcf population.unfiltered.vcf.gz \
  --cultivar-list cultivar.list \
  --wild-list wild.list \
  --keep-chr Chr01 \
  --outdir test_Chr01 \
  --prefix CA59_Chr01
```

Use a comma-separated list for multiple chromosomes, for example `--keep-chr Chr01,Chr02`.

## 6. Algorithm

### 6.1 One-pass genotype cache

On the first run, V17 reads the VCF once and stores SNP positions, `int8` genotypes, cultivar allele counts, and cultivar MAF by chromosome. The default cache is:

```text
OUTDIR/01_genotype_cache/
```

`OUTDIR/00_run_info/sample_qc.tsv` reports call rate and heterozygosity for every selected sample. The cache manifest records the VCF path, size, modification time, and sample lists. `--resume` validates these values before reusing the cache.

### 6.2 Common-SNP backbone

The default backbone views use cultivar MAF thresholds of 0.05 and 0.10. Common variants emphasize older, stable population structure, and agreement between the two views provides extra boundary support.

### 6.3 Multiple spatial scales

Default coarse scans are 1 Mb with a 500 kb step and 500 kb with a 250 kb step. Default fine scans are 100, 50, 25, and 10 kb. The 10 and 25 kb scans use 100 kb of context by default and expand to 200 kb when SNP information is insufficient. This provides 10 kb boundary resolution without treating every sparse 10 kb segment as a biological breakpoint.

### 6.4 Distance and clustering

The primary missing-aware distance is:

```text
distance_snp_per_mb = mismatch_count / compared_site_count
                      × eligible_site_count / block_length_Mb
```

The default backbone threshold is 100 SNP/Mb. Complete-linkage clustering is cut at this threshold, so the farthest comparable pair in every formal major group must pass it. The pipeline verifies the maximum within-group distance again after clustering.

With `--genotype-model inbred`, heterozygous calls are treated as missing; this is appropriate for highly homozygous lines. With `--genotype-model dosage`, 0/1/2 genotype states are retained for unphased diploids. For heterozygous populations, use `dosage` and inspect PCA and callability carefully.

### 6.5 Parent blocks, short gaps, and child blocks

Adjacent windows merge only when sample membership is stable; the default ARI requirement is 0.85. A coarse low-information gap is bounded by both:

- `--coarse-max-gap-bp`, default 1 Mb; and
- `--max-gap-fraction`, default 5% of the candidate block length.

Bridged regions are listed in `*.internal_gaps.tsv`. A sustained membership change creates a boundary. Fine scans then run inside every parent; child blocks are retained without deleting their parent block.

### 6.6 Whole-block validation

Windows propose boundaries. Every formal block is reclustered and checked across its full interval. By default, a coarse candidate must receive support from at least two window scales. This prevents an isolated local call from becoming a formal block.

### 6.7 Recent subhaplotypes

Within each major backbone haplotype, V17 recalculates group-specific MAF and uses variants polymorphic in that group to detect recent structure. The default recent threshold is 20 SNP/Mb, stricter than the backbone threshold.

- sufficiently frequent groups are named `SubHapXX`;
- smaller groups are retained as `MinorSubHapXX`;
- groups with too little SNP information are `UnresolvedSubHap` and are excluded from recent-subhap association candidates.

### 6.8 PCA

PCA is a QC and visualization layer. Missing values are mean-imputed in the PCA copy, variants are standardized, and points are colored by the final complete-linkage assignment. PCA does not overwrite the formal haplotype calls.

## 7. Important defaults

| Option | Default | Meaning |
|---|---:|---|
| `--backbone-mafs` | `0.05,0.10` | cultivar MAF thresholds for backbone calls |
| `--coarse-windows` | `1000000,500000` | parent-block windows |
| `--fine-windows` | `100000,50000,25000,10000` | child/boundary windows |
| `--backbone-threshold-snp-per-mb` | `100` | maximum within-backbone distance |
| `--recent-threshold-snp-per-mb` | `20` | maximum within-recent-subhap distance |
| `--membership-ari` | `0.85` | membership agreement across windows |
| `--min-scale-support` | `2` | minimum coarse window-scale support |
| `--min-structural-hap-frequency` | `0.05` | minimum frequency for a structural major haplotype |
| `--min-structural-hap-samples` | `3` | minimum samples for a structural major haplotype |
| `--min-association-hap-frequency` | `0.10` | minimum association haplotype frequency |
| `--min-association-assigned-fraction` | `0.50` | minimum combined major-group coverage |
| `--min-association-haplotypes` | `2` | minimum major groups in an association-ready block |
| `--max-window-snps` | `2000` | deterministic SNP cap per window distance calculation |
| `--max-block-snps` | `10000` | deterministic SNP cap during whole-block validation |

When a cap is exceeded, V17 samples evenly across genomic position. The operation is deterministic for identical inputs and parameters.

## 8. Parameter guidance

### 8.1 Very long centromeric blocks

Start with the default coarse settings. If a multi-megabase block is repeatedly broken by short low-information intervals, try a moderate increase:

```bash
--coarse-max-gap-bp 2000000 --max-gap-fraction 0.08
```

Inspect `*.internal_gaps.tsv`, `LOCAL_PASS_FRACTION`, PCA, and the membership on both sides before increasing the gap further. The fractional limit prevents a short candidate from crossing a disproportionately large interval.

### 8.2 Emphasizing 10 kb resolution

The default run already includes 10 kb scans. With high SNP density, context can be reduced:

```bash
--fine-context-bp 50000 --fine-max-context-bp 100000
```

When SNP density is low, retain the default 100/200 kb context instead of reducing `--min-window-snps` to a very small value.

### 8.3 Inbred material

Use the default `--genotype-model inbred`. Review `sample_qc.tsv`; unexpectedly high heterozygosity should trigger checks of sample identity, mapping, and genotype calling before distance thresholds are relaxed.

### 8.4 Heterozygous diploid material

Use `--genotype-model dosage`. This compares 0/1/2 states but does not separate the two haplotypes in a `0|1` genotype. True chromosome-level analysis requires prior phasing and a dedicated phased-haplotype implementation; V17.0.0 does not provide that model.

### 8.5 Association criteria

The requested default rule is implemented directly: each major haplotype is at least 10%, major groups jointly cover at least 50% of cultivars, and at least two major groups exist.

For a stricter 15%/70% filter:

```bash
--min-association-hap-frequency 0.15 \
--min-association-assigned-fraction 0.70
```

These settings control entry into association matrices. Smaller haplotypes remain in the complete frequency and sample tables.

## 9. Output files

### 9.1 Blocks and hierarchy

- `PREFIX.blocks.tsv`: main table containing both `COARSE` and `FINE` records.
- `PREFIX.blocks.bed`: BED intervals.
- `PREFIX.block_hierarchy.tsv`: fine-to-parent relationships.
- `PREFIX.internal_gaps.tsv`: bridged low-information or discordant intervals.
- `PREFIX.window_calls.tsv`: detailed status for every scale, MAF view, and window.

`START0` is zero-based and `END` is end-exclusive, following BED convention. Human-readable coordinates run from `START0 + 1` through `END`.

Useful columns in `PREFIX.blocks.tsv` include:

- `N_MAJOR_HAPLOTYPES`: groups passing the structural 5%/sample-count rule;
- `ASSIGNED_FRACTION`: cultivar fraction assigned to those structural major groups;
- `MAX_WITHIN_DISTANCE_SNP_PER_MB`: farthest pair within a major group;
- `MIN_BETWEEN_DISTANCE_SNP_PER_MB`: nearest pair in different major groups;
- `BETWEEN_WITHIN_RATIO`: minimum between-group distance divided by maximum within-group distance;
- `PAIR_CALLABLE_FRACTION`: proportion of sample pairs with comparable sites;
- `LOCAL_PASS_FRACTION`: proportion of constituent local calls that passed;
- `ASSOCIATION_READY`: whether the 10%/50% association criteria pass.

### 9.2 Haplotypes and carriers

- `PREFIX.block_haplotype_frequency.tsv`
- `PREFIX.block_haplotype_samples.tsv`
- `PREFIX.sample_by_block_haplotype.tsv`

`HapXX` identifies a structural major group and `MinorXX` identifies a smaller group. All samples remain in complete tables.

### 9.3 Recent subhaplotypes

- `PREFIX.backbone_to_recent_subhap.tsv`
- `PREFIX.recent_subhaplotype_frequency.tsv`
- `PREFIX.sample_by_recent_subhaplotype.tsv`
- `PREFIX.recent_subhap_association_status.tsv`
- `PREFIX.association_recent_subhap_onehot_matrix.tsv`

### 9.4 Association-ready outputs

- `PREFIX.association_ready_blocks.tsv`
- `PREFIX.association_excluded_blocks.tsv`, including failure reasons;
- `PREFIX.association_block_factor_matrix.tsv`, with one categorical variable per block;
- `PREFIX.association_haplotype_onehot_matrix.tsv`, with one binary column per major haplotype.

A block-level factor test is usually a useful first analysis, followed by haplotype-specific contrasts for significant blocks. Association models should also account for population structure, kinship, and experimental environment. V17 supplies predictors; it does not replace a complete GWAS or mixed-model framework.

### 9.5 PCA, figures, and reports

- `PREFIX.block_pca_scores.tsv`
- `PREFIX.block_pca_variance.tsv`
- `PREFIX.block_pca_top_loadings.tsv`
- `PREFIX.block_pca_plots.pdf`
- `PREFIX.block_landscape.pdf/png`
- `PREFIX.report.md`
- `PREFIX.output_index.tsv`

### 9.6 Wild projection

These files are generated only with `--wild-list`:

- `PREFIX.wild_to_cultivar_haplotype.tsv`: distances from each wild sample to every cultivar haplotype and the closest group;
- `PREFIX.wild_support_summary.tsv`: the number of wild samples closest to each cultivar haplotype.

Wild proximity is evidence of ancestral compatibility, not proof of direct ancestry. V17.0.0 does not automatically polarize alleles as ancestral or derived.

## 10. Using the output with phenotypes

Join a phenotype table to `PREFIX.association_block_factor_matrix.tsv` by `SAMPLE`. Each block is a categorical predictor. Use the one-hot matrix when the downstream tool requires numeric predictors.

Before modeling, verify:

1. the intersection of phenotype and VCF samples;
2. the observed carrier count of every tested haplotype;
3. whether haplotypes correlate strongly with genome-wide population structure;
4. whether coarse or fine blocks should be tested to control the multiple-testing burden;
5. whether a parent and its children are redundant before putting all of them in one ordinary linear model.

## 11. Performance and resource use

The dominant computation grows approximately with the square of the sample count. The main genotype cache requires roughly:

```text
retained SNP count × selected sample count × 1 byte
```

plus coordinate and frequency arrays. Recommended optimization order:

1. tune parameters on one chromosome using `--keep-chr`;
2. let the first run finish the cache, then use `--resume` for the same VCF and lists;
3. use `--skip-pca` and/or `--skip-recent-subhaps` for a quick block-only exploration;
4. for many samples, explore with `--max-window-snps 1000 --max-block-snps 5000`;
5. validate final results with larger SNP caps;
6. avoid running many large jobs against the same storage simultaneously.

Reuse a validated cache with:

```bash
python3 hapblock_v17.py [the original options] --resume
```

If the VCF or lists changed, cache validation stops the run. A new output directory is safest. To intentionally replace an output/cache, use `--force`; this deletes the specified existing output/cache, so check the path first.

## 12. Assessing accuracy

Review these items before biological interpretation:

1. `00_run_info/sample_qc.tsv` for low call rate or unusually high heterozygosity;
2. `blocks.tsv` for within-threshold groups and clear between-group separation;
3. `window_calls.tsv` for sustained boundaries rather than isolated low-information windows;
4. `internal_gaps.tsv` for blocks crossing excessively large anomalous intervals;
5. PCA plots for separation, gradients, or single-sample outliers;
6. haplotype frequency tables for adequate carrier counts;
7. parent/child redundancy before association testing.

After the first local run, the most useful files for a joint accuracy review are:

- `PREFIX.report.md`
- `PREFIX.blocks.tsv`
- `PREFIX.block_haplotype_frequency.tsv`
- `PREFIX.internal_gaps.tsv`
- `PREFIX.window_calls.tsv`
- `00_run_info/sample_qc.tsv`
- `PREFIX.block_landscape.png`
- both wild projection tables, when applicable.

## 13. Troubleshooting

### `Requested samples absent from VCF`

One or more list IDs do not exactly match the VCF header. Check case, hidden characters, and sample renaming.

### `VCF is not chromosome-sorted`

A chromosome occurs in multiple separated runs or positions decrease within a chromosome. Sort the VCF before running.

### Too few blocks

Check common-SNP counts, sample call rates, `--min-window-snps`, the distance threshold, and `--min-scale-support`. Inspect window failure states on one chromosome before substantially relaxing thresholds.

### Blocks are too fragmented

Determine whether low SNP density causes the breaks. Fine context or gap limits can be increased moderately. A sustained membership change should remain a true boundary.

### Blocks are too long

Inspect `internal_gaps.tsv` and `LOCAL_PASS_FRACTION`. Reduce `--coarse-max-gap-bp` or `--max-gap-fraction`, raise `--membership-ari`, and compare membership on both sides.

### Memory use is too high

Analyze by chromosome or reduce the selected samples. Lower SNP caps reduce computation and some temporary memory, but pairwise sample matrices still grow quadratically with sample count.

## 14. Tests and acknowledgment

Run the bundled synthetic tests with:

```bash
bash tests/run_tests.sh
```

The tests cover chain-merging prevention, short-gap bridging, sustained boundaries, the 10%/50% association rule, cultivar-only mode, and cultivar+wild mode. They do not run the user's large biological dataset.

OpenAI Codex assisted with requirement synthesis, workflow redesign, implementation, testing, and documentation for HapBlock V17. Biological conclusions should be evaluated by the researchers using species genetics, sampling design, and independent evidence.
