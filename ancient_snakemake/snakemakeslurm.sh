#! /bin/bash 
## defaults, will be overwritten if redefined in options parsing
script_name=$(basename $0)
prog_version=0.1.0

# default behavior: run all steps:
onlycmt=false
onlyposteager=false
dryrun=

## functions

# Function to display script usage
usage() {
    echo "Usage: $0 -I <file> -O <file> [options]"
    echo "Options:"
    echo "  -c, --cmt           Run only cmt snakemake (case)"
    echo "  -p, --posteager     Run only post eager snakemake (mapping)"
    echo "  -d, --dryrun        Dry run of snakemake (posteager/mapping)"
    echo "  -v, --version       Show version"
    echo "  -h, --help          Show help"
    exit 1
}
# Parse command-line options
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -c|--cmt)
            onlycmt=true
            shift
            ;;
        -p|--posteager)
            onlyposteager=true
            shift
            ;;
        -d|--dryrun)
            dryrun=--dry-run
            shift
            ;;
        -v|--version)
            version=true
            shift
            ;;
        -h|--help)
            help=true
            shift
            ;;
        *)
            echo "Unknown option: $key"
            usage
            ;;
    esac
done

help() { # print help, explanation for all parameters
    printf  "
    $(basename ${script_name^^})
    $script_name - Initialize unified ancient snakemake (post eager --> CMT generation)
    
    OPTIONS
        Options:
        -p
            --posteager - Only run the posteager.Snakefile. Only without also invoking -c/--cmt
        -c
            --cmt - Only run the cmt.Snakefile, if restarting a run that failed on this section. Only without also invoking -p/--posteager
        -d, 
            --dryrun - Dry run of snakemake (posteager/mapping)
        -h     
            --help - Print this help message
        -v      
            --version - Print version

    AUTHOR
        Ian Light-Maka (ilight1542@gmail.com)

    VERSION
        ${prog_version}
    \n    
    "
}

## validation ##
if [[ $version == true ]]; then
  printf $prog_version
  exit 0
fi

if [[ $help == true ]]; then help && exit 0; fi

if [[ $onlycmt == true && $onlyposteager == true  ]]
    then
    echo "ERROR: Only -c OR -p can be set at once"
    usage
fi


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

call_snakemake() {

  snakemake -p \
    --snakefile "$1" \
    --latency-wait 60 \
    -j 100 \
    --cluster-config "$(dirname "$0")/cluster.slurm.json" \
    --cluster "sbatch ${SM_ARGS}" \
    --cluster-status scripts/slurm_status.py \
    --default-resources "tmpdir='/ptmp/${USER}/tmp'" \
    --rerun-incomplete \
    --restart-times 1 \
    --keep-going \
    --use-conda \
    --conda-prefix /nexus/posix0/MPIIB-keylab/snakemake_conda_envs/ \
    --group-components make_link_group=100000 var2pos=200 \
    "${@:2}"

}

if [[ $onlyposteager == true ]]
then
    call_snakemake posteager.Snakefile ${dryrun}
    exit 0
fi

if [[ $onlycmt == true ]]
then
    call_snakemake cmt.Snakefile ${dryrun}
    exit 0
fi

# unified call (default behavior)
call_snakemake posteager.Snakefile && call_snakemake cmt.Snakefile
exit 0
