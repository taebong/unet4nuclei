experiment_name="03"   #options: 03 (for U2OS);01_MEF (MEF, or other nuclei with speckles)
input_dir="/n/groups/mitchison/Tae/SpinningDisk_IF/20210813_cTY36_ALFA-ORF6_mutants/After_fixation_100x/"
output_dir="/n/groups/mitchison/Tae/analysis/unet4nuclei_outputs/20210813_cTY36_ALFA-ORF6_mutants/After_fixation_100x/"
image_list_pth=$output_dir"seg_image_list.csv"
binning_factor=10
#pattern=".*_w.*642_.*.TIF"
pattern=".*405.tif"
#pattern=".*_DAPI.tif"

python generateImageLists_general.py $input_dir $output_dir $pattern
python batch_preprocessing.py $image_list_pth $input_dir $output_dir $binning_factor
python batch_prediction.py $experiment_name $image_list_pth $output_dir
