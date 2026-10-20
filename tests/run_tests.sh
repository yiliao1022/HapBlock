#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
python3 "$root/tests/test_cluster.py"
python3 "$root/tests/test_block_logic.py"
python3 "$root/tests/make_synthetic_v17.py"
python3 "$root/tests/test_cache.py"
bash "$root/tests/test_end_to_end.sh"
bash "$root/tests/test_cultivar_only.sh"
python3 -m py_compile "$root/hapblock_v17.py" "$root"/hapblock_v17lib/*.py
echo "PASS all HapBlock V17 tests"
