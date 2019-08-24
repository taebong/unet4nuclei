#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-12:00:00
#SBATCH -p short
#SBATCH --mem=8000
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
input_dir="/n/groups/mitchison/Tae/SpinningDisk/20190811_MEF_transport_glyco_perturbation/Well8/"
output_dir="/n/groups/mitchison/Tae/analysis/unet4nuclei_outputs/20190811_MEF_transport_glyco_perturbation/Well8/"

python nuclei_LEXY_analysis.py $input_dir $output_dir
python LEXY_model_fitting.py $output_dir

#done
