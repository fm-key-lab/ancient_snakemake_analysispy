# Local analysis initial QC and minimal analysis script

## Refactoring of code base to reduce the overwhelm and messiness
## additional goals:
- ease rerunning a similar but slightly adjusted filtering parameters, without creating a whole new script
- reporducing the same analysis
- moving filtering parameters out of script-by-script for easy comparison across runs
- standardize output directories

run within either an interactive session or a slurm script on logan (or raven, but logan is more appropriate IMO)

`
python /path/to/local_analysis_initial_qc/local_analysis_initial_qc_main.py -p test_local.json
`