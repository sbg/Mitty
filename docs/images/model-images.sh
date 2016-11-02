#!/usr/bin/env bash
set -x
for MODEL in "1kg-pcr-free" "hiseq-2500-v1-pcr-free" "hiseq-X-v1-HLI" "hiseq-X-v2.5-Garvan" "old-Garvan"
do
  mitty describe-read-model ${MODEL}.pkl ${MODEL}.png
done