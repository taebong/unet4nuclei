rootpth="/n/groups/mitchison/Tae/SpinningDisk/20191223_cTY10siOGA_cTY52siOGA_fusion/"
#rootpth="/n/groups/mitchison/Tae/WideField_NuclearTransport/20190321_intensity_test_20x/"
input_dir_pattern="Well*/"

pthlists=$(ls -d $rootpth$input_dir_pattern)

output_base="/n/groups/mitchison/Tae/analysis/unet4nuclei_outputs/20191223_cTY10siOGA_cTY52siOGA_fusion/"


#experiment_name="01_MEF"  #No need to change unless you want to change trained model
experiment_name="03"    #03 for U2OS, 01_MEF for MEF

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

python generateImageLists_wGFP.py $input_dir $output_dir

python batch_preprocessing.py $image_list_pth $input_dir $output_dir $binning_factor

python batch_prediction.py $experiment_name $image_list_pth $output_dir

python nuclei_LEXY_analysis_wGFP.py $input_dir $output_dir

python LEXY_model_fitting.py $output_dir

done
