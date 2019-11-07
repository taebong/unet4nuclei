import sys
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit

import pandas as pd
import os
import gc

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import colors as c

if len(sys.argv) != 2:
    print("Use: python LEXY_model_fitting.py output_dir")
    sys.exit()
else:
    analysis_dir = sys.argv[1]    #analysis pth

df_data = pd.read_csv(analysis_dir+'nuclei_tracking_results.csv')
#df_data.head(20)

transition_times = df_data.groupby('AcqState')['micro_T'].min().reset_index().set_index('AcqState')
if 'PostLit' not in transition_times.index:    #if PostLit images were not analyzed
    transition_times = transition_times.append(pd.DataFrame([(np.inf)],columns=['micro_T'],index=['PostLit']))

def LEXY_model(t,kout,kin,b,c,transition_times=transition_times):
    a = 1
    t1 = transition_times.at['Lit','micro_T']
    t2 = transition_times.at['PostLit','micro_T']
    y = (t<t1)*a   #PreLit
    y = y + ((t>=t1) & (t<t2))*((a-b)*np.exp(-(t-t1)*kout)+b)   #Lit
    
    bottom = a*np.exp(-(t2-t1)*kout)+b
   
    if ~np.isinf(t2): 
        y = y + (t>=t2)*((c-bottom)*(1-np.exp(-(t-t2)*kin))+bottom)     #PostLit
                     
    return y


def LINuS_model(t,kout,kin,b,c,transition_times=transition_times):
    a = 1
    t1 = transition_times.at['Lit','micro_T']
    t2 = transition_times.at['PostLit','micro_T']
    y = (t<t1)*a   #PreLit
    y = y + ((t>=t1) & (t<t2))*((b-a)*(1-np.exp(-(t-t1)*kin))+a)   #Lit
    
    top = (b-a)*(1-np.exp(-(t2-t1)*kin))+a
    
    if ~np.isinf(t2):
        y = y + (t>=t2)*(c+(top-c)*np.exp(-(t-t2)*kout))     #PostLit
   
    return y


def fit_model(df,model=None,
                   p0=None,
                  bounds=None,
                  log_file=None):
    grp = df.groupby(['Pos','Cycle','ID'])
   
    Ngrps = len(grp)
    cols = ['Pos','ID','Cycle','kout','kin','b','c','delta_kout','delta_kin','delta_b','delta_c','ChiSq','converged','N_dp']
    df_res = pd.DataFrame(columns=cols,index=np.arange(Ngrps))
    for i,g in enumerate(grp):
        ind = g[0]
        df_select = g[1]
        #remove na and zeros
        df_select = df_select[['micro_T','meanint_LEXY_normed','stdint_LEXY_normed','area']].replace([np.inf, -np.inf],np.nan).dropna()
        df_select[['meanint_LEXY_normed','stdint_LEXY_normed','area']] = df_select[['meanint_LEXY_normed','stdint_LEXY_normed','area']].replace([0],np.nan).dropna()  
        
        df_res.at[i,'Pos'] = ind[0]
        df_res.at[i,'Cycle'] = ind[1]
        df_res.at[i,'ID'] = ind[2]
        df_res.at[i,'N_dp'] = df_select.shape[0]
        
        if df_select.shape[0]==0:
            continue
        
        x = df_select['micro_T']
        y = df_select['meanint_LEXY_normed']
        dy = df_select['stdint_LEXY_normed']/np.sqrt(df_select['area'])  #sem
       
        #popt,pcov = curve_fit(model,x,y,p0=p0,sigma=dy,bounds=bounds) 
        try:
            popt,pcov = curve_fit(model,x,y,p0=p0,sigma=dy,bounds=bounds)
            
            dof = len(x)-len(popt) #degrees of freedom
            chisq = ((y-model(x,*popt))**2/dy**2).sum()/dof   #reduced chi_sq
        
            df_res.at[i,'converged'] = True
            df_res.at[i,['kout','kin','b','c']] = popt
            df_res.at[i,['delta_kout','delta_kin','delta_b','delta_c']] = np.sqrt(pcov.diagonal())  #errors
            df_res.at[i,'ChiSq'] = chisq

        except RuntimeError:   #Fit not converged
            df_res.at[i,'converged'] = False
            
        except:
            message = 'error raised for Pos %d, Cycle %d, ID %d' %ind
            print(message)
            if log_file:
                log_file.write(message+"\n")
            
    
    df_res = df_res.set_index(['Pos','Cycle','ID'])
    df_res[['kout','kin','b','c','delta_kout','delta_kin','delta_b','delta_c','ChiSq']]=df_res[['kout','kin','b','c','delta_kout','delta_kin','delta_b','delta_c','ChiSq']].apply(pd.to_numeric)
    
    return df_res


# Determine whether it is LEXY or LINuS
grp = df_data.groupby('AcqState')
means = grp['meanint_LEXY_normed'].mean()

if (means['Lit']<means['PreLit']):
    model = LEXY_model
    p0 = [1/10,1/30,0.7,1]
    bounds = ([0,0,0,0.5],[np.inf,np.inf,1,np.inf])
else:
    model = LINuS_model
    p0 = [1/10,1/30,2,1]
    bounds = ([0,0,1,0.5],[np.inf,np.inf,np.inf,np.inf])

print("Fitting %s....." %model.__name__)

log_file = open(analysis_dir+"fitting_logfile.txt",'w') 
log_file.write("Fitting model: "+model.__name__+"\n")

df_res = fit_model(df_data,model=model,p0=p0,bounds=bounds,log_file=log_file)

# merge with prelit average
df_prelit = df_data.loc[df_data['AcqState']=='PreLit']
grp = df_prelit.groupby(['Pos','Cycle','ID'])
var_list = ['eccen','area','perimeter','meanint_nucl','meanint_LEXY']
prelit_avg = grp[var_list].mean().reset_index().set_index(['Pos','Cycle','ID'])
df_res = df_res.join(prelit_avg)
#df_res.head()

# merge with other channel average (e.g. 447 avearage)
var_list = [v for v in df_data.columns if ('meanint_' in v) & ('nucl' not in v) & ('LEXY' not in v)]
grp = df_data.groupby(['Pos','Cycle','ID'])
avg = grp[var_list].mean().reset_index().set_index(['Pos','Cycle','ID'])
df_res = df_res.join(avg)
# save results
df_res.to_csv(analysis_dir+'model_fit_results.csv')

log_file.close()

print('Saving data and fitted model plots....')
os.makedirs(analysis_dir+'model_fit_plots/',exist_ok=True)
savedir  = analysis_dir+'model_fit_plots/'

#check the results
grp = df_res.reset_index().groupby(['Pos','Cycle'])
fig = plt.figure(figsize=(8,7));
curves_per_fig = 10

t1 = transition_times.at['Lit','micro_T']
t2 = transition_times.at['PostLit','micro_T']

for ind,df in grp:
    pos = ind[0]
    cyc = ind[1]
    df = df.dropna()

    count = 1
    fignum = 1
    for i,row_res in df.iterrows():
        nucl_id = row_res['ID']
        popt = tuple(row_res[['kout','kin','b','c']].values)

        df_select = df_data.loc[(df_data['Pos']==pos) & (df_data['Cycle']==cyc) & (df_data['ID']==nucl_id)]
        
        x = df_select['micro_T']
        y = df_select['meanint_LEXY_normed']
        dy = df_select['stdint_LEXY_normed']/np.sqrt(df_select['area'])  #sem

        plt.scatter(x,y,alpha=0.6,label='Data, %d' %nucl_id);

        xx = np.linspace(x.min(),x.max(),num=100)
        yy = model(xx,*popt)

        if row_res['converged']:
            color = 'k'
        else:
            color = 'r'

        plt.plot(xx,yy,color,label='Fit, %d' %nucl_id);

        if (count==curves_per_fig):
            plt.legend(title='Nucleus ID',ncol=2)
            plt.axvline(x=t1,linestyle='--',color='gray',alpha=0.6)
            
            if ~np.isinf(t2):
                plt.axvline(x=t2,linestyle='--',color='gray',alpha=0.6)
            
            plt.xlabel('Time (sec)')
            plt.ylabel('Norm. Intensity (A.U.)')

            fig.savefig(savedir+'Pos%02d_Cycle%03d_%03d.jpg' %(pos,cyc,fignum),
               frameon=False,facecolor=None,edgecolor=None);
            count=1
            fignum+=1
            
            plt.clf()
            gc.collect()
        else:
            count+=1


