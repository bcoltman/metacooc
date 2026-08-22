#!/usr/bin/env bash
set -euo pipefail

METACOOC_BIN="${METACOOC_BIN:-metacooc}"
DATA_RELEASE="R226_globdb_rev1"
OUTPUT_DIR="${1:-results/examples/bee-gut-association}"
data_args=(--data_release "$DATA_RELEASE")
if [[ -n "${METACOOC_DATA_DIR:-}" ]]; then
    data_args+=(--data_dir "$METACOOC_DATA_DIR")
fi

"$METACOOC_BIN" download "${data_args[@]}" --include_metadata
"$METACOOC_BIN" association \
    "${data_args[@]}" \
    --search_mode metadata \
    --search_string "bee gut" \
    --output_dir "$OUTPUT_DIR" \
    --filter_rank species \
    --min_conditional_probability 0.0 \
    --null_model FE \
    --q_threshold 0.10 \
    --label_top_n 6
