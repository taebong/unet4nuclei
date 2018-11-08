#!/bin/bash
source setup_modules.sh

experiment_name="20181029_LongTL"
input_dir="/n/groups/mitchison/Tae/SpinningDisk/20181029_LEXY_20xDry_LongTL/LongTL/"
output_dir="/home/ty118/analysis/unet4nuclei_outputs/20181029_LongTL/"
image_list_pth="/home/ty118/analysis/unet4nuclei_outputs/20181029_LongTL/seg_image_list.csv"

python batch_prediction.py $experiment_name $image_list_pth $input_dir $output_dir
