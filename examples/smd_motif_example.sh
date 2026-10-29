#!/usr/bin/env bash
set -euo pipefail

# Demonstration only. This short motif is not a full-length native protein model.
slime-proteo-md smd \
  --sequence PSEPGYPGQPSEPGYP \
  --conformation extended \
  --steps 2000 \
  --report-interval 20 \
  --pull-velocity 0.02 \
  --spring-k 500 \
  --temperature 300 \
  --platform CPU \
  --out results/smd_demo/
