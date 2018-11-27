rootpth="/n/groups/mitchison/Tae/SpinningDisk/20181116_cTY48_20xDry/"
input_dir_pattern="*"
output_base="/home/ty118/analysis/unet4nuclei_outputs/20181116_"

experiment_name="02_modified_model"  #No need to change unless you want to change trained model
binning_factor=2
listname="seg_image_list.csv"

for pth in $(ls -d $rootpth$input_dir_pattern)
do
input_dir=$pth"/" 
output_dir=$output_base$(basename $pth)"/"
mkdir -p $output_base$(basename $pth)

#check directory names
echo $input_dir
echo $output_dir

image_list_pth="$output_dir$listname"

python generateImageLists.py $input_dir $output_dir

python batch_preprocessing.py $image_list_pth $input_dir $output_dir $binning_factor

python batch_prediction.py $experiment_name $image_list_pth $output_dir

python nuclei_LEXY_analysis.py $input_dir $output_dir

python LEXY_model_fitting.py $output_dir

done
