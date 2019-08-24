experiment_name="03"
input_dir="/n/groups/mitchison/Tae/LongTermTL/20190711_cTY65_glyco_perturbation/converted_tiffs/"
output_dir="/n/groups/mitchison/Tae/analysis/unet4nuclei_outputs/20190711_cTY65_glyco_perturbation/"
#output_dir="/home/ty118/analysis/OMX_photoconversion/20190510_U2OS_photoconversion_glyco_perturbation/"
image_list_pth=$output_dir"seg_image_list.csv"
binning_factor=2
pattern=".*Far-Red_.*TIF"

python generateImageLists_general.py $input_dir $output_dir $pattern
python batch_preprocessing.py $image_list_pth $input_dir $output_dir $binning_factor
python batch_prediction.py $experiment_name $image_list_pth $output_dir
