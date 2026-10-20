#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
testdir="$root/tests"
python3 "$testdir/make_synthetic_v17.py"
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT
export MPLCONFIGDIR="$tmpdir/matplotlib"
export XDG_CACHE_HOME="$tmpdir/cache"
python3 "$root/hapblock_v17.py" \
  --vcf "$testdir/synthetic.vcf" \
  --cultivar-list "$testdir/cultivar.list" \
  --wild-list "$testdir/wild.list" \
  --outdir "$tmpdir/out" \
  --prefix test \
  --genotype-model inbred \
  --backbone-mafs 0.05,0.10 \
  --coarse-windows 200000,100000 \
  --fine-windows 50000,25000,10000 \
  --fine-context-bp 50000 \
  --fine-max-context-bp 100000 \
  --backbone-threshold-snp-per-mb 20 \
  --recent-threshold-snp-per-mb 10 \
  --min-window-snps 3 \
  --min-coarse-block-bp 100000 \
  --min-fine-block-bp 10000 \
  --coarse-max-gap-bp 100000 \
  --fine-max-gap-bp 20000 \
  --min-scale-support 2 \
  --pca-max-blocks 5
test -s "$tmpdir/out/test.blocks.tsv"
test -s "$tmpdir/out/test.block_haplotype_samples.tsv"
test -s "$tmpdir/out/test.association_ready_blocks.tsv"
test -s "$tmpdir/out/test.wild_to_cultivar_haplotype.tsv"
test -s "$tmpdir/out/test.report.md"
test -s "$tmpdir/out/00_run_info/sample_qc.tsv"
python3 - "$tmpdir/out" <<'PY'
import csv, pathlib, sys
out=pathlib.Path(sys.argv[1])
with (out/'test.blocks.tsv').open() as f: blocks=list(csv.DictReader(f,delimiter='\t'))
assert blocks, 'no blocks'
assert any(x['ASSOCIATION_READY']=='YES' for x in blocks), blocks
with (out/'test.block_haplotype_frequency.tsv').open() as f: freq=list(csv.DictReader(f,delimiter='\t'))
assert any(float(x['FREQUENCY_ALL_CULTIVARS']) >= .10 for x in freq)
with (out/'test.output_index.tsv').open() as f: index=list(csv.DictReader(f,delimiter='\t'))
assert any(x['OUTPUT']=='WILD_PROJECTION' and x['STATUS']=='PRESENT' for x in index)
print('PASS end-to-end',len(blocks),'blocks')
PY
