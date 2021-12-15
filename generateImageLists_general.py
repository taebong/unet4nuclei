import sys
import numpy as np
import pandas as pd

import os
import re

nsysarg = len(sys.argv)
print(sys.argv)
if (nsysarg != 3) & (nsysarg != 4):
    print("Use: python generateImageLists_general.py input_dir output_dir pattern(optional)")
    sys.exit()
else:
    pth = sys.argv[1]    #raw image pth
    savedir = sys.argv[2]   #analysis pth
    if nsysarg==4:
        pattern = sys.argv[3]
        regexp = re.compile(pattern)
    else:
        #pattern = '.*.tif'
        pattern = '(PreLit|Lit|PostLitGFP|PostLit|PostLitLSSmKate)\d*_w.*(642|491|447)\s{0,1}_.*.TIF'
        regexp = re.compile(pattern,re.IGNORECASE) #inclde all tif files

#print(sys.argv[3])
#print(pattern)
os.makedirs(savedir,exist_ok=True)

allfnames = os.listdir(pth)
fnames = [regexp.match(s).group(0) for s in allfnames if regexp.match(s)]

df = pd.DataFrame(fnames,columns=['DNA'])
#print(df)

#save image list
df.to_csv('%sseg_image_list.csv' %savedir,index=False)
