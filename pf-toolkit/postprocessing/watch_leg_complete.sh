#!/usr/bin/env bash
# Wait for a MicroSim leg's final frame to be COMPLETELY written, then run the
# interface-amplitude diagnostic over the whole leg.
#
# Why the size check matters: a frame that is still being written is short, and
# every VTK reader will happily parse the truncated file and give you a plausible
# but wrong answer. This waits for the exact expected byte count, twice, 60 s
# apart, before trusting it.
#
# Expected size for a MicroSim restart/frame file is
#     5 * 8 * nx * ny * nz + len(header) + 4
# (phi, mu, c blocks plus temperature; see reembed_frame.py for the layout).
#
# Usage, detached on a login node:
#   nohup bash watch_leg_complete.sh \
#       --jobid 12345 --frame /path/to/DATA/Processor_0/case_1500000.vtk \
#       --expected-bytes 283852816 --case-dir /path/to/case \
#       --outdir analysis_interface > watch.log 2>&1 &
#
# Optional:
#   --conda-sh /path/to/conda.sh --conda-env microsim-pp
#   --dx 1.0771e-8 --dt 2.0e-9 --substrate-top 20
#   --timeout-hours 9   --poll-seconds 120
set -u

JOBID=""; FRAME=""; EXPECTED_BYTES=""; CASE_DIR="."; OUTDIR="analysis_interface"
CONDA_SH=""; CONDA_ENV=""
DX=1.0771e-8; DT=2.0e-9; SUBSTRATE_TOP=20
TIMEOUT_HOURS=9; POLL=120

while [[ $# -gt 0 ]]; do
    case "$1" in
        --jobid)          JOBID=$2; shift 2 ;;
        --frame)          FRAME=$2; shift 2 ;;
        --expected-bytes) EXPECTED_BYTES=$2; shift 2 ;;
        --case-dir)       CASE_DIR=$2; shift 2 ;;
        --outdir)         OUTDIR=$2; shift 2 ;;
        --conda-sh)       CONDA_SH=$2; shift 2 ;;
        --conda-env)      CONDA_ENV=$2; shift 2 ;;
        --dx)             DX=$2; shift 2 ;;
        --dt)             DT=$2; shift 2 ;;
        --substrate-top)  SUBSTRATE_TOP=$2; shift 2 ;;
        --timeout-hours)  TIMEOUT_HOURS=$2; shift 2 ;;
        --poll-seconds)   POLL=$2; shift 2 ;;
        -h|--help)        sed -n '2,28p' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

: "${FRAME:?--frame is required}"
: "${EXPECTED_BYTES:?--expected-bytes is required}"

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DEADLINE=$(( $(date +%s) + TIMEOUT_HOURS * 3600 ))

log() { printf '%s %s\n' "$(date '+%F %T')" "$*"; }

file_size() { stat -c %s "$1" 2>/dev/null || stat -f %z "$1" 2>/dev/null; }

log "watching for $(basename "$FRAME") to reach ${EXPECTED_BYTES} bytes"

while :; do
    if [ "$(date +%s)" -gt "$DEADLINE" ]; then
        log "DEADLINE exceeded without a complete final frame - giving up"
        exit 2
    fi
    if [ -f "$FRAME" ] && [ "$(file_size "$FRAME")" = "$EXPECTED_BYTES" ]; then
        sleep 60
        if [ "$(file_size "$FRAME")" = "$EXPECTED_BYTES" ]; then
            log "final frame complete (${EXPECTED_BYTES} bytes)"
            break
        fi
        log "size was transient - still writing, continuing to wait"
    fi
    sleep "$POLL"
done

if [ -n "$JOBID" ]; then
    state=$(squeue -j "$JOBID" -h -o %T 2>/dev/null || true)
    log "job $JOBID state: ${state:-absent from queue}"
fi

if [ -n "$CONDA_SH" ]; then
    # shellcheck disable=SC1090
    source "$CONDA_SH"
    [ -n "$CONDA_ENV" ] && conda activate "$CONDA_ENV"
fi

cd "$CASE_DIR" || exit 1

# The solver creates shift.dat only once the moving window first fires. Before
# that the diagnostic must be told so explicitly; after it, absolute heights are
# wrong without it. RMS / peak-to-valley / spectrum are shift-invariant either way.
if [ -f DATA/shift.dat ]; then
    SHIFT_ARGS=(--shift-file DATA/shift.dat)
    log "shift.dat present - window has fired, using lab-frame heights"
else
    SHIFT_ARGS=(--no-shift)
    log "no shift.dat - asserting the window has not fired"
fi

log "running interface_amplitude.py over the full leg"
python "$HERE/interface_amplitude.py" DATA/Processor_0/*.vtk \
    --dx "$DX" --dt "$DT" --substrate-top "$SUBSTRATE_TOP" \
    "${SHIFT_ARGS[@]}" \
    --output-dir "$OUTDIR"
rc=$?
log "diagnostic exit $rc"
exit $rc
