#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-3:00:00
#SBATCH -p short
#SBATCH --mem=4000
#SBATCH -o log/lexyanalysis_%j.out
#SBTACH -e log/lexyanalysis_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=taeyeon_yoo@hms.harvard.edu


module load gcc/6.2.0
module load cuda/9.0
module load python/3.6.0

source nucleiseg_env/bin/activate


#for i in 1 2 3
#do
# SPECIFY THIS!
input_dir="/n/groups/mitchison/Tae/SpinningDisk/20190322_cTY48_hyperosmotic_shock/Control_Well1/"
output_dir="/home/ty118/analysis/unet4nuclei_outputs/20190329_widefield_analysis_addition_test/Control_Well1/"

python nuclei_LEXY_analysis.py $input_dir $output_dir

#done
