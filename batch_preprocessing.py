import sys
import os
import os.path

import numpy as np
import pandas as pd

from tqdm import tqdm

import skimage.io


if len(sys.argv) != 5:
    print("Use: python batch_preprocessing.py image_list.csv input_dir output_dir BINNING_FACTOR")
    sys.exit()
else:
    image_list = pd.read_csv(sys.argv[1])
    input_dir = sys.argv[2]
    output_dir = sys.argv[3]
    BINNING_FACTOR = int(sys.argv[4])

# Create output directories for transformed data

os.makedirs(output_dir+"normalized_images", exist_ok=True)

# # Image Preprocessing

filelist = image_list['DNA']

# run over all raw images
for filename in tqdm(filelist):

    # load image and its annotation
    orig_img = skimage.io.imread(input_dir + filename)       

    # binning
    sz = np.array(orig_img.shape)

    # if sz[0] and sz[1] are not multiples of BINNING_FACTOR, reduce them to the largest multiple of BINNING_FACTOR and crop image
    newsz = (sz/BINNING_FACTOR).astype(int)
    cropsz = newsz*BINNING_FACTOR
    orig_img = orig_img[0:cropsz[0],0:cropsz[1]]

    orig_img = orig_img.reshape((newsz[0],BINNING_FACTOR,newsz[1],BINNING_FACTOR))
    orig_img = orig_img.mean(-1).mean(1)

    # normalize to [0,1]
    percentile = 99.9
    high = np.percentile(orig_img, percentile)
    low = np.percentile(orig_img, 100-percentile)

    img = np.minimum(high, orig_img)
    img = np.maximum(low, img)

    img = (img - low) / (high - low) # gives float64, thus cast to 8 bit later
    img = skimage.img_as_ubyte(img) 

    skimage.io.imsave(output_dir+"normalized_images/" + filename[:-3] + 'png', img)    


print(img.dtype, img.shape)

