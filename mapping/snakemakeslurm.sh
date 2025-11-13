#!/bin/bash

# launch snakemake to run jobs via SLURM
# ntasks
SM_PARAMS="job-name ntasks time mem mail-user mail-type error output"

#SM_ARGS=" --parsable --cpus-per-task {cluster.cpus-per-task} --mem-per-cpu {cluster.mem-per-cpu-mb}"
SM_ARGS=" --parsable --cpus-per-task {cluster.cpus-per-task}" #new

for P in ${SM_PARAMS}; do SM_ARGS="${SM_ARGS} --$P {cluster.$P}"; done
echo "SM_ARGS: ${SM_ARGS}"

# our SLURM error/output paths expect a logs/ subdir in PWD
mkdir -p logs

conda config --set ssl_verify no

## activate snakemake conda env
source /u/iclight/bin/miniconda3/bin/activate snakemake

## Check if slurm_status.py is executable and change if not
if [[ ! -x "scripts/slurm_status.py" ]]
then
    chmod +x scripts/slurm_status.py;
    echo "Changed 'scripts/slurm_status.py' to executable";
fi

snakemake -p \
    $* \
    --latency-wait 60 \
    -j 100 \
    --cluster-config $(dirname $0)/cluster.slurm.json \
    --cluster "sbatch ${SM_ARGS}" \
    --cluster-status scripts/slurm_status.py \
    --default-resources "tmpdir='/ptmp/iclight/tmp'" \
    --rerun-incomplete \
    --restart-times 3 \
    --keep-going \
    --use-conda \
    --conda-prefix /nexus/posix0/MPIIB-keylab/snakemake_conda_envs/ \
    --group-components make_link_group=100000

    #--dry-run \
    #--group-components \
    # --dag \
    # | dot -Tsvg > dag.svg
