experiment_name="03"
output_dir="/n/groups/mitchison/Tae/analysis/unet4nuclei_outputs/Patrick_seg_test/"
image_list_pth=$output_dir"seg_image_list.csv"

python batch_prediction.py $experiment_name $image_list_pth $output_dir
