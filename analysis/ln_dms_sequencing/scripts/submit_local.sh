snakemake \
    --jobs 999 \
    --use-singularity \
    --singularity-args "--bind /fh --bind /hpc" \
    --cores=all \
    $@
