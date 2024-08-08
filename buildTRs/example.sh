#!/bin/bash

snakemake -s buildTRs.smk \
    --configfile example.yaml \
    --cluster "{params.cluster} -o {params.stdout} -e {params.stderr}" \
    --jobname "{params.rule}.{jobid}.sh" \
    --jobs 200 \
    -d run \
    -p \
    --restart-times 1 \
    --keep-going \
    --rerun-incomplete --dry-run
