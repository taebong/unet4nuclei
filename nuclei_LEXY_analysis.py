from utils_LEXY import *

import sys
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit

import pandas as pd
import skimage
from skimage import io
from skimage import measure
import mahotas as mh

import matplotlib as mpl
from matplotlib import pyplot as plt

import trackpy as tp

import os,glob
import re
import gc
import json

if len(sys.argv) != 3:
    print("Use: python nuclei_LEXY_analysis.py raw_image_dir analysis_dir")
    sys.exit()
else:
    raw_image_dir = sys.argv[1]    #raw image pth
    analysis_dir = sys.argv[2]   #analysis pth

code_ver = 'v5'
refinement_setting = {'N':3, 'repeat': 1}   #segmentation refinement setting. 
#N=3, repeat=1 for rupture assay, and N=10, repeat=2 for import/export assay 
tracking_setting = {'link_distance':40,'memory':1,'adaptive_step':0.99,'adaptive_stop':5}     #nucleus tracking setting
    
basename = ''
AcqStates = ['PreLitZScan','PreLit','Lit','PostLit','PostLitGFP','Rupture','PostLitLSSmKate']
#reg_cycle = '(?<=[PreLitZScan|PreLit|Lit|PostLit|Rupture])(?P<Cycle>\d*)(?=_w)'
reg_Pos = '(?<=_s)(?P<Pos>\d*)'
reg_T = '(?<=_t)(?P<T>\d*)(?=.)'
reg_Ch = '(?<=_w\d)(tae |Fluo |Confocal )?(?P<Ch>[^\s^_]*)\s?\w+?(?=_)'
drop_Chs = ['Cyan']
#nucl_Chs = ['642','Far-Red']
#lexy_Chs = ['561','Red']
Ch_map = {'642':'nucl',
          'Far-Red':'nucl',
          '640':'nucl',
          '561':'LEXY',
          'Red':'LEXY',
          '447':'447',
          '491':'491',
          '445':'445',
          'LSSmKate':'LSSmKate'}  # Ch:ChannelName dictionary. Channel names should include 'nucl' at least  
#binning_factor = 2

#Load data list

fpths = glob.glob(raw_image_dir+basename+'*.TIF')+glob.glob(raw_image_dir+basename+'*.tif')

metadict_list = []
for fpth in fpths:
    fname = re.split('/',fpth)[-1]
    
    metadata = {'Filename':fname,'AcqState':np.nan,'T':0,'Ch':np.nan,'Cycle':0,'Pos':0}
    
    state_re = re.search('(?<=%s).*?(?=\d*_)' %basename,fname)
    T_re = re.search(reg_T,fname)
    Ch_re = re.search(reg_Ch,fname)
    cyc_re = re.search('(?P<Cycle>\d*)(?=_)',fname)
    Pos_re = re.search(reg_Pos,fname)
    
    if state_re:
        metadata['AcqState'] = state_re.group(0)
    
    if T_re:
        metadata['T'] = int(T_re.group('T'))
        
    if Ch_re:
        metadata['Ch'] = Ch_re.group('Ch')
        
    if cyc_re:
        if len(cyc_re.group('Cycle'))>0:
            metadata['Cycle'] = int(cyc_re.group('Cycle'))
            
    if Pos_re:
        if len(Pos_re.group('Pos'))>0:
            metadata['Pos'] = int(Pos_re.group('Pos'))
        
    metadict_list.append(metadata)
    
df_meta = pd.DataFrame(metadict_list)
df_meta['AcqState'] = pd.Categorical(df_meta['AcqState'],AcqStates)
df_meta[['T','Cycle','Pos']] = df_meta[['T','Cycle','Pos']].astype(int,errors='ignore')
df_meta = df_meta.loc[[ch not in drop_Chs for ch in df_meta['Ch']]].copy()
df_meta = df_meta.sort_values(['Pos','Cycle','AcqState','T','Ch']).reset_index(drop=True)

df_meta.head()


#Load segmentation list
fpths = glob.glob(analysis_dir+'segm/'+basename+'*.png')

dict_list = []
for fpth in fpths:
    fname = re.split('/',fpth)[-1]
    
    d = {'Filename':fname,'AcqState':np.nan,'T':0,'Ch':np.nan,'Cycle':0,'Pos':0}
    
    state_re = re.search('(?<=%s).*?(?=\d*_)' %basename,fname)
    T_re = re.search(reg_T,fname)
    Ch_re = re.search(reg_Ch,fname)
    cyc_re = re.search('(?<=%s)(?P<Cycle>\d*)(?=_)' %basename,fname)
    Pos_re = re.search(reg_Pos,fname)
    
    if state_re:
        d['AcqState'] = state_re.group(0)
    
    if T_re:
        d['T'] = int(T_re.group('T'))
        
    if Ch_re:
        d['Ch'] = Ch_re.group('Ch')
        
    if cyc_re:
        if len(cyc_re.group('Cycle'))>0:
            d['Cycle'] = int(cyc_re.group('Cycle'))
            
    if Pos_re:
        if len(Pos_re.group('Pos'))>0:
            d['Pos'] = int(Pos_re.group('Pos'))
        
    dict_list.append(d)
    
df_seg = pd.DataFrame(dict_list)
df_seg['AcqState'] = pd.Categorical(df_seg['AcqState'],AcqStates)
df_seg[['T','Cycle','Pos']] = df_seg[['T','Cycle','Pos']].astype(int,errors='ignore')
df_seg = df_seg.loc[[ch not in drop_Chs for ch in df_seg['Ch']]].copy()
df_seg = df_seg.sort_values(['Pos','Cycle','AcqState','T','Ch']).reset_index(drop=True)

df_seg.head()


# refine segmentation: take intersection of N consecutive images, and watershed. Do this forward and backward and repeat
print('refining segmentation....')
os.makedirs(analysis_dir+'segm_refined/',exist_ok=True)

N = refinement_setting['N']
repeat = refinement_setting['repeat']

same_previous_setting = False
if os.path.isfile(analysis_dir+'segm_refined/refinement_setting.json'):    
    try:
        with open(analysis_dir+'segm_refined/refinement_setting.json','r') as fp:
            m = json.load(fp)
        same_previous_setting = (m['Code_ver']==code_ver) & (m['N']==N) & (m['repeat']==repeat)
    except:
        print("previous refinement setting cannot be loaded. Rerunning refinement..")

if same_previous_setting:
    print('refined segmentation already exists')
else:
    df_grp = df_seg.groupby(['Pos','Cycle'])
    for ind,df in df_grp:
        print("Pos %d, Cycle %d processing" %ind)

        #load segmentation results
        labels = []

        for i,row in df.reset_index().iterrows():
            label = io.imread(analysis_dir+'segm/'+row['Filename'])
            labels.append(label)
        labels = np.array(labels)
        segs = labels>0

        #Foward and backward
        k = 0
        while (k<2*(repeat+1)):
            for i in range(df.shape[0]-N):
                intersect = np.prod(segs[i:i+N],axis=0)
                intersect_l = measure.label(intersect)
                distances = mh.distance(segs[i])
                surface = (distances.max() - distances)

                areas = mh.cwatershed(surface,intersect_l)
                new_label = areas*segs[i]

                labels[i] = new_label
                segs[i] = new_label>0

            labels = labels[::-1]
            segs = segs[::-1]
            k += 1

        #Save refined segmentations
        for i,row in df.reset_index().iterrows():
            io.imsave(analysis_dir+'segm_refined/'+row['Filename'],labels[i])

    print('Done! \n')

    # save refining setup
    m = {'Code_ver':code_ver,'N':N,'repeat':repeat}

    with open(analysis_dir+'segm_refined/refinement_setting.json','w') as fp:
        json.dump(m,fp)

        

print("extracting features from raw images and segmentation results....")
df_data = pd.DataFrame()
lexy_bg = 1e5 #background level
Ch_map_filtered = {ch:Ch_map[ch] for ch in Ch_map if ch in np.unique(df_meta['Ch'].values)}

current_pos = -99
for i,row in df_seg.iterrows():
    if current_pos != row['Pos']:
        fr = 0    #reset frame number to 0 for each field
        current_pos = row['Pos']
    
    label=io.imread(analysis_dir+'segm_refined/'+row['Filename'])

    print(row['Filename'])
    
    props_dict = {}
    macro_T_dict = {}
    for select_ch in Ch_map_filtered:
        ch_name = Ch_map_filtered[select_ch]
        
        ind = (df_meta['Cycle']==row['Cycle']) & (df_meta['Pos']==row['Pos']) & (df_meta['AcqState']==row['AcqState']) & (df_meta['Ch']==select_ch) & (df_meta['T']==row['T'])
        
        if ind.sum()==0:
            continue
        else:
            fname = df_meta.at[np.where(ind)[0][0],'Filename']
        
        img = io.imread(os.path.join(raw_image_dir,fname))
        
        if 'binning_factor' not in locals():
            binning_factor = int(img.shape[0]/label.shape[0])
        
        img = downsample(img,binning_factor)
        
        if (ch_name == 'LEXY') & (img.min() <lexy_bg):
            lexy_bg = img.min()
            
        props_dict[ch_name] = measure.regionprops(label,intensity_image=img,coordinates='rc')
        
        macro_T_dict[ch_name] = get_time(os.path.join(raw_image_dir,fname))
    
    newlist = []

    if 'nucl' not in props_dict:
        print('nucl channel not exist for AcqState %s, Cycle %d, Pos %d, T %d' %(row['AcqState'],row['Cycle'],row['Pos'],row['T']))
        continue
    
    for j,p_nucl in enumerate(props_dict['nucl']):
        #nucleus dimensions
        cent = p_nucl.centroid
        orientation = p_nucl.orientation
        area = p_nucl.area
        eccen = p_nucl.eccentricity
        perimeter = p_nucl.perimeter
        
        macro_T = macro_T_dict['nucl']

        newdict = {'AcqState':row['AcqState'],'Cycle':row['Cycle'],
                   'Pos':row['Pos'],'pretrack_ID':p_nucl.label,
                   'cent_x':cent[1],'cent_y':cent[0],
                   'orientation':orientation,'area':area,
                   'perimeter':perimeter,'eccen':eccen,
                   'T':row['T'],'macro_T':macro_T,'frame':fr}
        
        # add channel-specific properties to newdict
        for ch_name in props_dict:
            props = props_dict[ch_name]
            p = props[j]
            
            intim = p.intensity_image
            
            newdict['meanint_'+ch_name] = p.mean_intensity
            newdict['maxint_'+ch_name] = p.max_intensity
            newdict['minint_'+ch_name] = p.min_intensity
            newdict['stdint_'+ch_name] = np.std(intim[intim>0])
        
        newlist.append(newdict)
        
    df_data = pd.concat((df_data,pd.DataFrame(newlist)))
    
    fr += 1
    
    del img, label

df_data = df_data.reset_index(drop=True)

print("Done!\n")


#Track nuclei
pos_list = df_data['Pos'].unique()

link_distance = tracking_setting['link_distance']
memory = tracking_setting['memory']
adaptive_step = tracking_setting['adaptive_step']
adaptive_stop = tracking_setting['adaptive_stop']

df_data_tracked = pd.DataFrame()
for pos in pos_list:
    print("tracking nuclei for Pos %d...." %int(pos))
    
    df_select = df_data.loc[df_data['Pos']==pos].copy().reset_index(drop=True)
    
    df_select = tp.link(df_select,link_distance,
                        pos_columns=['cent_x','cent_y'],
                        t_column='frame',memory=memory,
                        adaptive_step=adaptive_step,adaptive_stop=adaptive_stop)
    df_select = df_select.rename(columns={'particle':'ID'})
    df_select = df_select.sort_values(['frame','ID'])
    df_select['ID'] = df_select['ID']+1   # to make it start from 1

    df_data_tracked = pd.concat((df_data_tracked,df_select))
    
    print("Done!\n")

df_data = df_data_tracked.reset_index(drop=True)

del df_data_tracked
del df_select

# save tracking setup
m = {'Code_ver':code_ver,'link_distance':link_distance,'memory':memory,'adaptive_step':adaptive_step,'adaptive_stop':adaptive_stop}

with open(analysis_dir+'tracking_setting.json','w') as fp:
    json.dump(m,fp)
    

#normalized meanint to PreLit average and set macro time (clock starts from the beginning of data set) 
# and micro time (clock starts from the beginning of each cycle)
# For rupture assay, no normalization

print("normalizing....")

grp = df_data.loc[df_data['AcqState']=='PreLit'].groupby(['Pos','Cycle','ID'])
var_list = ['eccen','area','perimeter']
var_list = var_list + ['meanint_'+ch_name for ch_name in Ch_map_filtered.values()]
prelit_avg = grp[var_list].mean().reset_index().set_index(['Pos','Cycle','ID'])

for ind,row in prelit_avg.iterrows():
    mask = (df_data['Pos']==ind[0]) & (df_data['Cycle']==ind[1]) & (df_data['ID']==ind[2])
    df_data.loc[mask,'meanint_LEXY_normed'] = (df_data.loc[mask,'meanint_LEXY']-lexy_bg)/(row['meanint_LEXY']-lexy_bg)
    df_data.loc[mask,'stdint_LEXY_normed'] = df_data.loc[mask,'stdint_LEXY']/(row['meanint_LEXY']-lexy_bg)
    #df_data.loc[mask,'meanint_LEXY_normed'] = (df_data.loc[mask,'meanint_LEXY']-lowest_val.at[i,'meanint_LEXY'])/(row['meanint_LEXY']-lowest_val.at[i,'meanint_LEXY'])

df_data['macro_T'] = df_data['macro_T'] - np.min(df_data['macro_T'])
    
grp_T = df_data.groupby(['Pos','Cycle'])['macro_T'].min().reset_index().set_index(['Pos','Cycle'])
for ind,row in grp_T.iterrows():
    mask = (df_data['Pos']==ind[0]) & (df_data['Cycle']==ind[1])
    df_data.loc[mask,'micro_T'] = df_data.loc[mask,'macro_T']-row['macro_T']
    

print("Done!\n")


# Save df_data and movies

print("Saving tracking data.....")
df_data.to_csv(analysis_dir+'nuclei_tracking_results.csv',index=False)

# savedir = analysis_dir+'nuclei_tracking_image/'
# os.makedirs(savedir,exist_ok=True)

# rmap = getLabelColorMap()

# fig = plt.figure();
# for i,row in df_seg.iterrows():
#     df_fr = df_data.loc[(df_data['Cycle']==row['Cycle']) 
#                         & (df_data['Pos']==row['Pos'])
#                         & (df_data['AcqState']==row['AcqState'])
#                        & (df_data['T']==row['T'])].copy()
    
#     label=io.imread(analysis_dir+'segm_refined/'+row['Filename'])
    
#     lexy_fname = df_meta.at[np.where((df_meta['Cycle']==row['Cycle'])
#                                      & (df_meta['Pos']==row['Pos'])
#                                      & (df_meta['AcqState']==row['AcqState']) 
#                                    & np.array([ch in lexy_Chs for ch in df_meta['Ch']]) 
#                                    & (df_meta['T']==row['T']))[0][0],'Filename']
#     lexy = io.imread(raw_image_dir+lexy_fname)
#     lexy = downsample(lexy,binning_factor)
    
#     nucl_fname = df_meta.at[np.where((df_meta['Cycle']==row['Cycle'])
#                                      & (df_meta['Pos']==row['Pos'])
#                                      & (df_meta['AcqState']==row['AcqState']) 
#                                    & np.array([ch in nucl_Chs for ch in df_meta['Ch']]) 
#                                    & (df_meta['T']==row['T']))[0][0],'Filename']
#     nucl = io.imread(raw_image_dir+nucl_fname)
#     nucl = downsample(nucl,binning_factor)
    
#     percentile = 99.9
#     if i == 0:  #normalization factors for images  
#         lexy_high = np.percentile(lexy,percentile)
#         lexy_low = np.percentile(lexy,100-percentile)
#         nucl_high = np.percentile(nucl,percentile)
#         nucl_low = np.percentile(nucl,100-percentile)
    
#     lexy = normalize_image(lexy,low=lexy_low,high=lexy_high)
#     nucl = normalize_image(nucl,low=nucl_low,high=nucl_high)
    
#     # change from pretrack label to posttrack label
#     id_map = df_fr[['pretrack_ID','ID']].drop_duplicates().set_index('pretrack_ID')
#     new_label = label.copy() #post track label

#     for preID in np.unique(label[label>0]):
#         new_label = np.where(label==preID,id_map.at[preID,'ID'],new_label)

#     fig,axs = showSegmentation(new_label,lexy,nucl,
#                                rmap,df_fr,fig=fig,
#                                t=df_fr['micro_T'].values[0],
#                               state=row['AcqState']);
    
#     fig.savefig(savedir+'frame%06d_Pos%02d_Cycle%03d_%s_%03d.jpg' %(i,int(row['Pos']),int(row['Cycle']),row['AcqState'],int(row['T'])),
#                frameon=False,facecolor=None,edgecolor=None,quality=80);

#     plt.clf()
    
#     del df_fr
#     del lexy, nucl, label
    
#     gc.collect()
    

print("Done!\n")


