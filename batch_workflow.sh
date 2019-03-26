rootpth="/n/groups/mitchison/Tae/SpinningDisk/20190322_cTY48_hyperosmotic_shock/"
input_dir_pattern="*/"

pthlists=$(ls -d $rootpth$input_dir_pattern)

output_base="/home/ty118/analysis/unet4nuclei_outputs/20190322_cTY48_hyperosmotic_shock/"

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

cp batch_workflow.sh $output_dir

python generateImageLists.py $input_dir $output_dir

python batch_preprocessing.py $image_list_pth $input_dir $output_dir $binning_factor

python batch_prediction.py $experiment_name $image_list_pth $output_dir

python nuclei_LEXY_analysis.py $input_dir $output_dir

python LEXY_model_fitting.py $output_dir

done
