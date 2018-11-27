input_dir="/n/groups/mitchison/Tae/SpinningDisk/20181029_LEXY_20xDry_LongTL/SingleTimePoint2/"
output_dir="/home/ty118/analysis/unet4nuclei_outputs/20181029_SingleTP2/"
image_list_pth="/home/ty118/analysis/unet4nuclei_outputs/20181029_SingleTP2/seg_image_list.csv"

python batch_preprocessing.py $image_list_pth $input_dir $output_dir 2
