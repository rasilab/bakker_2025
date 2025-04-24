snakemake \
    --jobs 999 \
    --latency-wait 60 \
    --cluster-config cluster.yaml \
    --cluster "sbatch -n {cluster.n}  -t {cluster.time}" \
    --use-conda \
    --use-singularity \
    --singularity-args "--bind /fh --bind /hpc" \
    --cores=all \
    -p $@
