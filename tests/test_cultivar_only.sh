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
  --outdir "$tmpdir/out" \
  --prefix cultivar_only \
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
  --skip-pca
test -s "$tmpdir/out/cultivar_only.blocks.tsv"
test -s "$tmpdir/out/00_run_info/sample_qc.tsv"
test ! -e "$tmpdir/out/cultivar_only.wild_to_cultivar_haplotype.tsv"
python3 - "$tmpdir/out" <<'PY'
import csv, pathlib, sys
out = pathlib.Path(sys.argv[1])
with (out / 'cultivar_only.output_index.tsv').open() as handle:
    index = {row['OUTPUT']: row['STATUS'] for row in csv.DictReader(handle, delimiter='\t')}
assert index['WILD_PROJECTION'] == 'NOT_GENERATED', index
report = (out / 'cultivar_only.report.md').read_text()
assert 'Mode: `CULTIVAR_ONLY`' in report
assert 'No ancestral/derived interpretation' in report
with (out / '00_run_info/sample_qc.tsv').open() as handle:
    qc = list(csv.DictReader(handle, delimiter='\t'))
assert len(qc) == 40 and all(row['GROUP'] == 'CULTIVAR' for row in qc)
print('PASS cultivar-only mode')
PY
