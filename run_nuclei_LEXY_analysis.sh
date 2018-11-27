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

for i in 1 2 3
do
# SPECIFY THIS!
input_dir="/n/groups/mitchison/Tae/SpinningDisk/20181029_LEXY_20xDry_LongTL/SingleTimePoint$i/"
output_dir="/home/ty118/analysis/unet4nuclei_outputs/20181029_SingleTP$i/"

python nuclei_LEXY_analysis.py $input_dir $output_dir

done
