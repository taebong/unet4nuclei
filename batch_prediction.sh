experiment_name="02_modified_model"
output_dir="/home/ty118/analysis/unet4nuclei_outputs/20181029_SingleTP2/"
image_list_pth="/home/ty118/analysis/unet4nuclei_outputs/20181029_SingleTP2/seg_image_list.csv"

python batch_prediction.py $experiment_name $image_list_pth $output_dir
