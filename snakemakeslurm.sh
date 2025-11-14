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
activate_conda_path=$(echo $CONDA_EXE | sed 's#bin/conda#bin/activate#g')
source ${activate_conda_path} snakemake

## Check if slurm_status.py is executable and change if not
if [[ ! -x "scripts/slurm_status_script.py" ]]
then
    chmod +x scripts/slurm_status_script.py;
    echo "Changed 'scripts/slurm_status_script.py' to executable";
fi

## Check if email has been updated
email_update=$(grep "UPDATE_EMAIL_ADDRESS" cluster.slurm.json)
[ ! ${#email_update} -eq 0 ] && echo "Please modify cluster.slurm.json to include your email address for error reporting" && exit 1

## create tmp and conda storage directory in /ptmp/ (if not already present)
mkdir -p /ptmp/${USER}/tmp
mkdir -p /nexus/posix0/MPIIB-keylab/snakemake_conda_envs/

## Check if slurm_status.py is executable and change if not
if [[ ! -x "scripts/slurm_status.py" ]]
then
    chmod +x scripts/slurm_status.py;
    echo "Changed 'scripts/slurm_status.py' to executable";
fi

snakemake -p \
    --snakefile posteager.Snakefile \
    $* \
    --latency-wait 60 \
    -j 100 \
    --cluster-config $(dirname $0)/cluster.slurm.json \
    --cluster "sbatch ${SM_ARGS}" \
    --cluster-status scripts/slurm_status.py \
    --default-resources "tmpdir='/ptmp/${USER}/tmp'" \
    --rerun-incomplete \
    --restart-times 1 \
    --keep-going \
    --use-conda \
    --conda-prefix /nexus/posix0/MPIIB-keylab/snakemake_conda_envs/ \
    --group-components make_link_group=100000 \
&& \
snakemake -p \
    --snakefile cmt.Snakefile \
    $* \
    --latency-wait 60 \
    -j 100 \
    --cluster-config $(dirname $0)/cluster.slurm.json \
    --cluster "sbatch ${SM_ARGS}" \
    --cluster-status scripts/slurm_status_script.py \
    --default-resources "tmpdir='/ptmp/${USER}/tmp'" \
    --group-components var2pos=200 \
    --rerun-incomplete \
    --restart-times 1 \
    --keep-going \
    --use-conda \
    --conda-prefix /nexus/posix0/MPIIB-keylab/snakemake_conda_envs/


    #--dry-run \
    #--group-components \
    # --dag \
    # | dot -Tsvg > dag.svg
