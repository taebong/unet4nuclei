rootpth="/n/groups/mitchison/Tae/SpinningDisk/20190723_MEF_transport_glyco_perturbation/"
input_dir_pattern="Well*/"

pthlists=$(ls -d $rootpth$input_dir_pattern)

output_base="/n/groups/mitchison/Tae/analysis/unet4nuclei_outputs/20190723_MEF_transport_glyco_perturbation/"

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

python generateImageLists.py $input_dir $output_dir

python batch_preprocessing.py $image_list_pth $input_dir $output_dir $binning_factor

done
