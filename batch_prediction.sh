experiment_name="03"
input_dir="/n/groups/mitchison/Tae/SpinningDisk/20190411_cTY48_glyco_perturbation_RL2/"
output_dir="/home/ty118/analysis/IF_analysis/20190411_cTY48_glyco_perturbation_RL2/"
image_list_pth=$output_dir"seg_image_list.csv"
binning_factor=2
patten=".*DAPI\.TIF"

python generateImageLists_general.py $input_dir $output_dir $pattern
python batch_preprocessing.py $image_list_pth $input_dir $output_dir $binning_factor
python batch_prediction.py $experiment_name $image_list_pth $output_dir
