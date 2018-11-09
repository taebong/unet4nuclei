#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-3:00:00
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=4000
#SBATCH -o log/batch_prediction_%j.out
#SBTACH -e log/batch_prediction_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=taeyeon_yoo@hms.harvard.edu

module load gcc/6.2.0
module load cuda/9.0
module load python/3.6.0

source nucleiseg_env/bin/activate

experiment_name="02_modified_model"
input_dir="/n/groups/mitchison/Tae/SpinningDisk/20181029_LEXY_20xDry_LongTL/LongTL/"
output_dir="/home/ty118/analysis/unet4nuclei_outputs/20181029_LongTL/"
image_list_pth="/home/ty118/analysis/unet4nuclei_outputs/20181029_LongTL/seg_image_list.csv"

python batch_prediction.py $experiment_name $image_list_pth $input_dir $output_dir
