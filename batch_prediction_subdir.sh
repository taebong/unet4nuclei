rootpth="/n/groups/mitchison/Tae/SpinningDisk_IF/20210804_GFP_concentration_measurement/"
input_dir_pattern="*/"
file_pattern=".*640.tif"

pthlists=$(ls -d $rootpth$input_dir_pattern)

output_base="/n/groups/mitchison/Tae/analysis/unet4nuclei_outputs/20210804_GFP_concentration_measurement_Station10/"

experiment_name="03"  #No need to change unless you want to change trained model
binning_factor=2
listname="seg_image_list.csv"

for pth in $pthlists
do
input_dir=$pth 
output_dir=$output_base$(basename $pth)"/"
mkdir -p $output_base$(basename $pth)

#check directory names
echo $input_dir
echo $output_dir

image_list_pth="$output_dir$listname"


python generateImageLists_general.py $input_dir $output_dir $file_pattern

python batch_preprocessing.py $image_list_pth $input_dir $output_dir $binning_factor

python batch_prediction.py $experiment_name $image_list_pth $output_dir


done
