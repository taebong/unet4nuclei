input_dir="/n/groups/mitchison/Tae/test_data/"
output_dir="/n/groups/mitchison/Tae/analysis/unet4nuclei_outputs/Patrick_seg_test/"
image_list_pth="/n/groups/mitchison/Tae/analysis/unet4nuclei_outputs/Patrick_seg_test/seg_image_list.csv"

python batch_preprocessing.py $image_list_pth $input_dir $output_dir 1
