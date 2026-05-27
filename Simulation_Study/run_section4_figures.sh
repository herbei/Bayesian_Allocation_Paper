#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INVOCATION_DIR="$(pwd)"
cd "$SCRIPT_DIR"

if [[ $# -ge 1 && -n "$1" ]]; then
  case "$1" in
    /*) CHAIN_FILE="$1" ;;
    *) CHAIN_FILE="$INVOCATION_DIR/$1" ;;
  esac
else
  CHAIN_FILE="$SCRIPT_DIR/Data/results/simulation_study_selected_run.RData"
fi

if [[ ! -f "$CHAIN_FILE" ]]; then
  echo "Missing MCMC artifact: $CHAIN_FILE" >&2
  echo "Run Rscript run_mcmc.R first, or pass an existing .RData artifact." >&2
  exit 1
fi

CHAIN_ABS="$(cd "$(dirname "$CHAIN_FILE")" && pwd)/$(basename "$CHAIN_FILE")"
FIGURE_ROOT="$SCRIPT_DIR/Figures"
OUTPUT_DIR="$FIGURE_ROOT/out/$(basename "${CHAIN_ABS%.RData}")"
mkdir -p "$OUTPUT_DIR"

FK_ENSEMBLE_FILE="$SCRIPT_DIR/Data/inputs/oxy_ensemble_20_members.rds" \
FK_OUTPUT_DIR="$FIGURE_ROOT/out" \
  Rscript "$SCRIPT_DIR/Figures/src/Figure01.R"

FK_ENSEMBLE_FILE="$SCRIPT_DIR/Data/inputs/oxy_ensemble_20_members.rds" \
FK_OUTPUT_DIR="$FIGURE_ROOT/out" \
  Rscript "$SCRIPT_DIR/Figures/src/Figure02.R"

FK_ENSEMBLE_FILE="$SCRIPT_DIR/Data/inputs/oxy_ensemble_20_members.rds" \
OXY_TRUE_FILE="$SCRIPT_DIR/Data/inputs/trueXX.txt" \
OXY_API_FILE="$SCRIPT_DIR/Data/inputs/api.txt" \
FK_OUTPUT_DIR="$FIGURE_ROOT/out" \
  Rscript "$SCRIPT_DIR/Figures/src/Figure03.R"

OXY_FIGURE_CHAIN_FILE="$CHAIN_ABS" \
OXY_FIGURE_OUTPUT_DIR="$OUTPUT_DIR" \
  Rscript "$SCRIPT_DIR/Figures/src/Figure04a.R"

OXY_FIGURE_CHAIN_FILE="$CHAIN_ABS" \
OXY_FIGURE_OUTPUT_DIR="$OUTPUT_DIR" \
  Rscript "$SCRIPT_DIR/Figures/src/Figure04b.R"

OXY_FIGURE_CHAIN_FILE="$CHAIN_ABS" \
OXY_FIGURE_OUTPUT_DIR="$OUTPUT_DIR" \
  Rscript "$SCRIPT_DIR/Figures/src/Figure05.R"

echo "Wrote static figures to $FIGURE_ROOT/out and posterior figures to $OUTPUT_DIR"
