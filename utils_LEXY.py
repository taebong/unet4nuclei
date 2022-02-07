import numpy as np

#from skimage.external import tifffile
import tifffile
from skimage import segmentation

import re
import datetime

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import colors as c


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
    try: 
        datestr = re.search('(?<=datetime \(\d{2}s\) b\').*(?=\.\d*\')',tf.info()).group(0)
        ts = datetime.datetime.strptime(datestr,'%Y%m%d %H:%M:%S').timestamp()
    except:
        datestr = tf.imagej_metadata['time']
        ts = datetime.datetime.strptime(datestr,'%Y-%m-%dT%H:%M:%S').timestamp()
    return ts



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

