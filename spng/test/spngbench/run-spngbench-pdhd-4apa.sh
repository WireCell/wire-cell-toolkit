#!/bin/bash
#
# Run the PDHD 4-APA spngbench grid and produce its cross-grid report.
#
# No arguments.  Everything is baked into spngbench-pdhd-4apa.json:
#   - napa = 4 (four per-APA pipelines; PDHD's physical maximum)
#   - device modes: CPU only, CPU + 1 GPU, and CPU + 2 GPU
#   - the 2-GPU case uses ONLY the "transverse" scheme, i.e. each APA pipeline
#     stays entirely on one GPU (no longitudinal sharding)
#
# Output goes to <work root>/spngbench-pdhd-4apa/ ; the grid report is written to
# <work root>/spngbench-pdhd-4apa/grid-index/.  The work root is the directory
# that contains the toolkit checkout (paths are derived, none are hard-wired).

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # .../toolkit/spng/test/spngbench
spngbench="$here/spngbench.py"
config="$here/spngbench-pdhd-4apa.json"

# Work root = parent of the toolkit checkout (here/../../../.. = <root>/toolkit's parent).
workroot="$(cd "$here/../../../.." && pwd)"
cd "$workroot"

outdir="spngbench-pdhd-4apa"

echo "== spngbench PDHD 4-APA grid =="
echo "   work root : $workroot"
echo "   config    : $config"
echo "   output    : $workroot/$outdir"
echo

# --mode all runs the cpu, single-GPU, and multi-GPU (shard) cell families.
python3 "$spngbench" grid --mode all --config "$config"

echo
echo "== grid report =="
python3 "$spngbench" grid-report "$outdir/grid-index.json" -o "$outdir/grid-index"

echo
echo "Done.  See $workroot/$outdir/grid-index/grid-summary.{md,html,tex}"
