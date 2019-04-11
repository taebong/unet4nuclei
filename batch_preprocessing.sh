input_dir="/n/groups/mitchison/Tae/SpinningDisk/20190404_cTY52_glyco_perturbation_RL2/"
output_dir="/home/ty118/analysis/IF_analysis/20190404_cTY52_glyco_perturbation_RL2/"
image_list_pth=$output_dir"seg_image_list.csv"
binning_factor=2

python batch_preprocessing.py $image_list_pth $input_dir $output_dir $binning_factor
