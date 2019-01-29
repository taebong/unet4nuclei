input_dir="/n/groups/mitchison/Tae/SpinningDisk/20190113_cTY48_cycle_arrest/Rupture/"

output_dir="/home/ty118/analysis/unet4nuclei_outputs/20190113_cTY48_Rupture/"

experiment_name="03"  #No need to change unless you want to change trained model
binning_factor=2
listname="seg_image_list.csv"

image_list_pth="$output_dir$listname"

python generateImageLists.py $input_dir $output_dir

python batch_preprocessing.py $image_list_pth $input_dir $output_dir $binning_factor

python batch_prediction.py $experiment_name $image_list_pth $output_dir

python nuclei_LEXY_analysis.py $input_dir $output_dir

#python rupture_analysis.py $input_dir $output_dir  #TBA


