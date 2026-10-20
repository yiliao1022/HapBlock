# HapBlock V17

HapBlock V17 detects hierarchical haplotype blocks from an unfiltered population VCF. Cultivar samples define all block boundaries and haplotype groups. An optional wild sample set is projected onto the completed cultivar result and never changes the cultivar calls.

Main features:

- common-SNP backbone blocks at multiple scales;
- coarse parent blocks plus fine children down to 10 kb boundary resolution;
- complete-linkage grouping with a verified maximum within-group SNP/Mb distance;
- recent subhaplotypes called within each common-SNP backbone haplotype;
- no artificial cap on the number of reported haplotypes;
- automatic cultivar-only and cultivar+wild modes;
- association-ready factor and one-hot matrices using the 10%/50% rule;
- one-pass binary genotype cache, PCA quality-control outputs, gap annotations, and run metadata.

## Install

```bash
conda env create -f environment.yml
conda activate hapblock-v17
```

or:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

## Run

Cultivar-only mode:

```bash
python3 hapblock_v17.py \
  --vcf population.unfiltered.vcf.gz \
  --cultivar-list cultivar.list \
  --outdir result_v17 \
  --prefix population
```

Cultivar+wild mode:

```bash
python3 hapblock_v17.py \
  --vcf population.unfiltered.vcf.gz \
  --cultivar-list cultivar.list \
  --wild-list wild.list \
  --outdir result_v17 \
  --prefix population
```

See [the Chinese manual](docs/HapBlock_V17_manual_CN.md) or [the English manual](docs/HapBlock_V17_manual_EN.md) before running a large dataset.

## Validate the installation

```bash
bash tests/run_tests.sh
```

The test suite uses generated synthetic data; it does not use the user's biological dataset.

## Version

This package is HapBlock V17.0.0. V16 files are not required and are not modified by this package.
