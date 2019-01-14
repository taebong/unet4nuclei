import sys
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit

import pandas as pd
import skimage
from skimage.external import tifffile
from skimage import segmentation
from skimage import io
from skimage import measure
import mahotas as mh

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import colors as c

import datetime

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

code_ver = 'v1.0'
refinement_setting = {'N':3, 'repeat': 1}   #segmentation refinement setting. 
#N=3, repeat=1 for rupture assay, and N=10, repeat=2 for import/export assay 
tracking_setting = {'link_distance':25,'memory':1}     #nucleus tracking setting
    
basename = ''
AcqStates = ['PreLitZScan','PreLit','Lit','PostLit','Rupture']
reg_cycle = '(?<=[PreLitZScan|PreLit|Lit|PostLit|Rupture])(?P<Cycle>\d*)(?=_w)'
reg_Pos = '(?<=_s)(?P<Pos>\d*)'
reg_T = '(?<=_t)(?P<T>\d*)(?=.)'
reg_Ch = '(?<=_w\d)(?P<Ch>\d{3})'
drop_Chs = ['447']
#binning_factor = 2

#Load data list

fpths = glob.glob(raw_image_dir+basename+'*.TIF')

metadict_list = []
for fpth in fpths:
    fname = re.split('/',fpth)[-1]
    
    metadata = {'Filename':fname,'AcqState':np.nan,'T':0,'Ch':np.nan,'Cycle':0,'Pos':0}
    
    state_re = re.search('(?<=%s).*?(?=\d*_)' %basename,fname)
    T_re = re.search(reg_T,fname)
    Ch_re = re.search(reg_Ch,fname)
    cyc_re = re.search(reg_cycle,fname)
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
    cyc_re = re.search(reg_cycle,fname)
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


def downsample(im,binning_factor):
    sz = np.array(im.shape)
    
    # if sz[0] and sz[1] are not multiples of BINNING_FACTOR, reduce them to the largest multiple of BINNING_FACTOR and crop image
    newsz = (sz/binning_factor).astype(int)
    cropsz = newsz*binning_factor
    im = im[0:cropsz[0],0:cropsz[1]]

    newim = im.reshape((newsz[0],binning_factor,newsz[1],binning_factor))
    return newim.mean(-1).mean(1)

def get_time(fpth):
    tf = tifffile.TiffFile(fpth)
    datestr = re.search('(?<=datetime \(\d{2}s\) b\').*(?=\.\d*\')',tf.info()).group(0)
    return datetime.datetime.strptime(datestr,'%Y%m%d %H:%M:%S').timestamp()




# refine segmentation: take intersection of N consecutive images, and watershed. Do this forward and backward and repeat
print('refining segmentation....')
os.makedirs(analysis_dir+'segm_refined/',exist_ok=True)

N = refinement_setting['N']
repeat = refinement_setting['repeat']

same_previous_setting = False
if os.path.isfile(analysis_dir+'segm_refined/refinement_setting.json'):
    with open(analysis_dir+'segm_refined/refinement_setting.json','r') as fp:
        m = json.load(fp)    
    same_previous_setting = (m['Code_ver']==code_ver) & (m['N']==N) & (m['repeat']==repeat)

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

current_pos = -99
for i,row in df_seg.iterrows():
    if current_pos != row['Pos']:
        fr = 0    #reset frame number to 0 for each field
        current_pos = row['Pos']
    
    label=io.imread(analysis_dir+'segm_refined/'+row['Filename'])
    
    lexy_fname = df_meta.at[np.where((df_meta['Cycle']==row['Cycle'])
                                     & (df_meta['Pos']==row['Pos'])
                                     & (df_meta['AcqState']==row['AcqState']) 
                                   & (df_meta['Ch']=='561') 
                                   & (df_meta['T']==row['T']))[0][0],'Filename']
    lexy = io.imread(raw_image_dir+lexy_fname)
    
    if 'binning_factor' not in locals():
        binning_factor = int(lexy.shape[0]/label.shape[0])

    lexy = downsample(lexy,binning_factor)
    
    if (lexy.min()<lexy_bg):
        lexy_bg = lexy.min()
    
    nucl_fname = df_meta.at[np.where((df_meta['Cycle']==row['Cycle']) 
                                     & (df_meta['Pos']==row['Pos'])
                                     & (df_meta['AcqState']==row['AcqState']) 
                                   & (df_meta['Ch']=='642') 
                                   & (df_meta['T']==row['T']))[0][0],'Filename']
    nucl = io.imread(raw_image_dir+nucl_fname)
    
    nucl = downsample(nucl,binning_factor)
    
    props_lexy = measure.regionprops(label,intensity_image=lexy,coordinates='rc')
    props_nucl = measure.regionprops(label,intensity_image=nucl,coordinates='rc')
    
    macro_T = get_time(raw_image_dir+lexy_fname)
    
    newlist = []
    for j in range(len(props_lexy)):
        p_lexy = props_lexy[j]
        p_nucl = props_nucl[j]
        
        cent = p_lexy.centroid
        orientation = p_lexy.orientation
        area = p_lexy.area
        eccen = p_lexy.eccentricity
        perimeter = p_lexy.perimeter
        meanint_LEXY = p_lexy.mean_intensity
        meanint_nucl = p_nucl.mean_intensity        
        maxint_LEXY = p_lexy.max_intensity        
        minint_LEXY = p_lexy.min_intensity  
    
        intim_LEXY = p_lexy.intensity_image
        intim_nucl = p_nucl.intensity_image
        
        stdint_LEXY = np.std(intim_LEXY[intim_LEXY>0])
        stdint_nucl = np.std(intim_nucl[intim_nucl>0])
        
        newlist.append({'AcqState':row['AcqState'],'Cycle':row['Cycle'],
                        'Pos':row['Pos'],
                        'pretrack_ID':p_nucl.label,
                         'cent_x':cent[1],'cent_y':cent[0],
                        'orientation':orientation,'area':area,
                        'perimeter':perimeter,
                        'eccen':eccen,'meanint_LEXY':meanint_LEXY,
                        'meanint_nucl':meanint_nucl,'maxint_LEXY':maxint_LEXY,
                        'minint_LEXY':minint_LEXY,
                        'stdint_LEXY':stdint_LEXY,'stdint_nucl':stdint_nucl,
                        'T':row['T'],'macro_T':macro_T,'frame':fr})
        
    df_data = pd.concat((df_data,pd.DataFrame(newlist)))
    
    fr += 1
    
    del nucl, lexy, label
    
df_data = df_data.reset_index(drop=True)

print("Done!\n")


#Track nuclei
pos_list = df_data['Pos'].unique()

link_distance = tracking_setting['link_distance']
memory = tracking_setting['memory']

df_data_tracked = pd.DataFrame()
for pos in pos_list:
    print("tracking nuclei for Pos %d...." %int(pos))
    
    df_select = df_data.loc[df_data['Pos']==pos].copy().reset_index(drop=True)
    
    df_select = tp.link(df_select,link_distance,
                        pos_columns=['cent_x','cent_y'],
                        t_column='frame',memory=memory)
    df_select = df_select.rename(columns={'particle':'ID'})
    df_select = df_select.sort_values(['frame','ID'])
    df_select['ID'] = df_select['ID']+1   # to make it start from 1

    df_data_tracked = pd.concat((df_data_tracked,df_select))
    
    print("Done!\n")

df_data = df_data_tracked.reset_index(drop=True)

del df_data_tracked
del df_select

# save tracking setup
m = {'Code_ver':code_ver,'link_distance':link_distance,'memory':memory}

with open(analysis_dir+'tracking_setting.json','w') as fp:
    json.dump(m,fp)
    

#normalized meanint to PreLit average and set macro time (clock starts from the beginning of data set) 
# and micro time (clock starts from the beginning of each cycle)
# For rupture assay, no normalization

print("normalizing....")

grp = df_data.loc[df_data['AcqState']=='PreLit'].groupby(['Pos','Cycle','ID'])
prelit_avg = grp[['meanint_LEXY','eccen','meanint_nucl','area','perimeter']].mean().reset_index().set_index(['Pos','Cycle','ID'])

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


def getLabelColorMap():
    colors = plt.cm.jet(range(256))
    np.random.shuffle(colors)
    colors[0] = (0.,0.,0.,1.)
    #rmap = c.ListedColormap(colors)
    return colors

def normalize_image(im,low=None,high=None):
    if low==None:
        low = np.min(im)
    if high==None:
        high = np.max(im)
    
    im = np.minimum(high,im)
    im = np.maximum(low,im)
    
    im = (im-low)/(high-low)
    
    return im    

def showSegmentation(label_im,norm_im1,norm_im2,rmap,df_track,zoom=3,fig=None,t=None,state=None):
    
    sz = label_im.shape
    
    combined = np.moveaxis([norm_im1,norm_im2,np.zeros(sz)],0,2)
    
    if fig:
        axs = fig.subplots(2, 2,sharex=True, sharey=True);
    else:
        fig,axs = plt.subplots(2, 2,sharex=True, sharey=True);
    
    w,h = plt.figaspect(sz[0]/sz[1]);
    fig.set_size_inches(w * zoom, h * zoom);
    
    axs[0][0].imshow(rmap[label_im%256]);    
    axs[0][1].imshow(segmentation.mark_boundaries(norm_im1,label_im,mode='inner',color=None,outline_color=[1,0,0]));
    axs[1][0].imshow(segmentation.mark_boundaries(norm_im2,label_im,mode='inner',color=None,outline_color=[1,0,0]));
    axs[1][1].imshow(segmentation.mark_boundaries(combined,label_im,mode='inner',color=None,outline_color=[1,1,1]));
    
    labels = np.unique(label_im)
    
    #X,Y = np.meshgrid(np.arange(sz[1]),np.arange(sz[0]))
    
    for l in labels[labels!=0]:
        #mask = (label_im == l)
        #xc = np.sum(X*mask)/np.sum(mask)
        #yc = np.sum(Y*mask)/np.sum(mask)
        
        xc,yc=df_track.loc[df_track['ID']==l,['cent_x','cent_y']].values[0]
        
        for ax in axs.flatten()[:3]:
            ax.text(xc,yc,l,fontsize=3*zoom,
                     horizontalalignment='center',
                     verticalalignment='center',color='k');
            
    for ax in axs.flatten():
        ax.set_yticks([]);
        ax.set_xticks([]);
    
    if t!=None: #add timepoint
        axs[0][0].text(sz[1]-10,10,'%d sec' %t,color='w',fontsize=7*zoom,horizontalalignment='right',
                       verticalalignment='top', bbox=dict(facecolor='black', alpha=0.5));
    
    if state: #add acquisition state
        axs[0][0].text(sz[1]-10,sz[0]-10,'%s' %state,color='w',fontsize=7*zoom,horizontalalignment='right',
                       verticalalignment='bottom', bbox=dict(facecolor='black', alpha=0.5));
    
    
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0.02,wspace=0.02);
    
    fig.tight_layout()

    return fig,axs



# Save df_data and movies

print("Saving tracking data.....")
df_data.to_csv(analysis_dir+'nuclei_tracking_results.csv',index=False)

savedir = analysis_dir+'nuclei_tracking_image/'
os.makedirs(savedir,exist_ok=True)

rmap = getLabelColorMap()

fig = plt.figure();
for i,row in df_seg.iterrows():
    df_fr = df_data.loc[(df_data['Cycle']==row['Cycle']) 
                        & (df_data['Pos']==row['Pos'])
                       & (df_data['T']==row['T'])].copy()
    
    label=io.imread(analysis_dir+'segm_refined/'+row['Filename'])
    
    lexy_fname = df_meta.at[np.where((df_meta['Cycle']==row['Cycle'])
                                     & (df_meta['Pos']==row['Pos'])
                                     & (df_meta['AcqState']==row['AcqState']) 
                                   & (df_meta['Ch']=='561') 
                                   & (df_meta['T']==row['T']))[0][0],'Filename']
    lexy = io.imread(raw_image_dir+lexy_fname)
    lexy = downsample(lexy,binning_factor)
    
    nucl_fname = df_meta.at[np.where((df_meta['Cycle']==row['Cycle'])
                                     & (df_meta['Pos']==row['Pos'])
                                     & (df_meta['AcqState']==row['AcqState']) 
                                   & (df_meta['Ch']=='642') 
                                   & (df_meta['T']==row['T']))[0][0],'Filename']
    nucl = io.imread(raw_image_dir+nucl_fname)
    nucl = downsample(nucl,binning_factor)
    
    percentile = 99.9
    if i == 0:  #normalization factors for images  
        lexy_high = np.percentile(lexy,percentile)
        lexy_low = np.percentile(lexy,100-percentile)
        nucl_high = np.percentile(nucl,percentile)
        nucl_low = np.percentile(nucl,100-percentile)
    
    lexy = normalize_image(lexy,low=lexy_low,high=lexy_high)
    nucl = normalize_image(nucl,low=nucl_low,high=nucl_high)
    
    # change from pretrack label to posttrack label
    id_map = df_fr[['pretrack_ID','ID']].drop_duplicates().set_index('pretrack_ID')
    new_label = label.copy() #post track label

    for preID in np.unique(label[label>0]):
        new_label = np.where(label==preID,id_map.at[preID,'ID'],new_label)

    fig,axs = showSegmentation(new_label,lexy,nucl,
                               rmap,df_fr,fig=fig,
                               t=df_fr['micro_T'].values[0],
                              state=row['AcqState']);
    
    fig.savefig(savedir+'frame%06d_Pos%02d_Cycle%03d_%s_%03d.jpg' %(i,int(row['Pos']),int(row['Cycle']),row['AcqState'],int(row['T'])),
               frameon=False,facecolor=None,edgecolor=None,quality=80);

    plt.clf()
    
    del df_fr
    del lexy, nucl, label
    
    gc.collect()
    

print("Done!\n")
