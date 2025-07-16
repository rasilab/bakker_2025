snakemake \
    --use-conda \
    --cores=8 \
    --use-singularity \
    --singularity-args "--bind /fh" \
    -p $@
