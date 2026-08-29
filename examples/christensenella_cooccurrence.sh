#!/usr/bin/env bash
set -euo pipefail

METACOOC_BIN="${METACOOC_BIN:-metacooc}"
DATA_RELEASE="R226_globdb_rev1"
OUTPUT_DIR="${1:-results/examples/smithii-cooccurrence}"
data_args=(--data_release "$DATA_RELEASE")
if [[ -n "${METACOOC_DATA_DIR:-}" ]]; then
    data_args+=(--data_dir "$METACOOC_DATA_DIR")
fi

"$METACOOC_BIN" download "${data_args[@]}" --include_metadata

common_args=(
    "${data_args[@]}"
    --search_mode focal_taxa
    --min_taxa_count 10
    --taxa_count_rank species
    --min_sample_count 50
    --null_model FE
    --q_threshold 0.10
    --null_scope metadata
    --null_metadata_query "human gut metagenome"
)

"$METACOOC_BIN" cooccurrence \
    "${common_args[@]}" \
    --search_string "s__Methanocatella smithii" \
    --filter_rank species \
    --min_conditional_probability 0.10 \
    --max_pairs 150000 \
    --label_top_n 5 \
    --output_dir "$OUTPUT_DIR/neighbourhood"

"$METACOOC_BIN" cooccurrence \
    "${common_args[@]}" \
    --search_string "g__Christensenella -> s__Methanocatella smithii" \
    --filter_rank species \
    --min_conditional_probability 0.0 \
    --max_pairs 5000000 \
    --label_top_n 3 \
    --output_dir "$OUTPUT_DIR/christensenella-species"

"$METACOOC_BIN" cooccurrence \
    "${common_args[@]}" \
    --aggregated \
    --search_string "g__Christensenella AGGREGATED -> s__Methanocatella smithii" \
    --min_conditional_probability 0.0 \
    --max_pairs 500000 \
    --label_top_n 1 \
    --output_dir "$OUTPUT_DIR/christensenella-aggregated"
