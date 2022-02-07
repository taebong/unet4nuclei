import sys
import numpy as np
import scipy as sp
import pandas as pd
import skimage
from skimage.external import tifffile

import os,glob
import re

if len(sys.argv) != 3:
    print("Use: python generateImageLists.py input_dir output_dir")
    sys.exit()
else:
    pth = sys.argv[1]    #raw image pth
    savedir = sys.argv[2]   #analysis pth

os.makedirs(savedir,exist_ok=True)

basename = ''
AcqStates = ['PreLitZScan','PreLit','Lit','PostLit','PostLitGFP','Rupture','PostLitLSSmKate']
reg_T = '(?<=_t)(?P<T>\d*)(?=.TIF)'
reg_Ch = '(?<=_w\d)(tae |Fluo |Confocal )?(?P<Ch>[^\s^_]*)\s?\w+?(?=_)'
drop_Chs = ['447','Cyan']

keep_Chs = ['642','Far-Red','640']
keep_AcqStates = ['PreLit','Lit','PostLit','PostLitGFP','Rupture','PostLitLSSmKate']
os.makedirs(savedir,exist_ok=True)


fpths = glob.glob(pth+basename+'*.TIF')

ImSeqs = dict()

metadict_list = []
for fpth in fpths:
    fname = re.split('/',fpth)[-1]
    
    metadata = {'Filename':fname,'AcqState':np.nan,'T':np.nan,'Ch':np.nan}
    
    state_re = re.search('(?<=%s).*?(?=\d*_)' %basename,fname)
    T_re = re.search(reg_T,fname)
    Ch_re = re.search(reg_Ch,fname)
    
    if state_re:
        metadata['AcqState'] = state_re.group(0)
    
    if T_re:
        metadata['T'] = int(T_re.group('T'))
        
    if Ch_re:
        metadata['Ch'] = Ch_re.group('Ch')
        
    metadict_list.append(metadata)
    
df_meta = pd.DataFrame(metadict_list)
df_meta['AcqState'] = pd.Categorical(df_meta['AcqState'],AcqStates)
df_meta = df_meta.loc[[ch not in drop_Chs for ch in df_meta['Ch']]].copy()
df_meta = df_meta.sort_values(['AcqState','T','Ch']).reset_index(drop=True)


df_select = df_meta.loc[[(c in keep_Chs) & (s in keep_AcqStates) for c,s in df_meta[['Ch','AcqState']].values]]
df_select = df_select.reset_index(drop=True)
df_select = df_select.rename(columns={'Filename':'DNA'})

#save image list
df_select.to_csv('%sseg_image_list.csv' %savedir,index=False)
