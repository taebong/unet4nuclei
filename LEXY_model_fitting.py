import sys
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit

import pandas as pd

import matplotlib as mpl
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

def LEXY_model(t,kout,kin,b,c,transition_times=transition_times):
    a = 1
    t1 = transition_times.at['Lit','micro_T']
    t2 = transition_times.at['PostLit','micro_T']
    y = (t<t1)*a   #PreLit
    y = y + ((t>=t1) & (t<t2))*((a-b)*np.exp(-(t-t1)*kout)+b)
    
    bottom = a*np.exp(-(t2-t1)*kout)+b
    
    y = y + (t>=t2)*((c-bottom)*(1-np.exp(-(t-t2)*kin))+bottom)
                     
    return y

def fit_LEXY_model(df,model=LEXY_model,
                   p0=[1/10,1/30,0.7,1],
                  bounds=([0,0,0,0.5],[np.inf,np.inf,1,np.inf])):
    grp = df.groupby(['Cycle','ID'])
   
    Ngrps = len(grp)
    cols = ['ID','Cycle','kout','kin','b','c','delta_kout','delta_kin','delta_b','delta_c','ChiSq','converged','N_dp']
    df_res = pd.DataFrame(columns=cols,index=np.arange(Ngrps))
    for i,g in enumerate(grp):
        ind = g[0]
        df_select = g[1]
        #remove na and zeros
        df_select = df_select[['micro_T','meanint_LEXY_normed','stdint_LEXY_normed','area']].replace([np.inf, -np.inf],np.nan).dropna()
        df_select[['meanint_LEXY_normed','stdint_LEXY_normed','area']] = df_select[['meanint_LEXY_normed','stdint_LEXY_normed','area']].replace([0],np.nan).dropna()  
        
        df_res.at[i,'Cycle'] = ind[0]
        df_res.at[i,'ID'] = ind[1]
        df_res.at[i,'N_dp'] = df_select.shape[0]
        
        if df_select.shape[0]==0:
            continue
        
        x = df_select['micro_T']
        y = df_select['meanint_LEXY_normed']
        dy = df_select['stdint_LEXY_normed']/np.sqrt(df_select['area'])  #sem
        
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
            print('error raised for Cycle %d, ID %d' %ind)
            
    
    df_res = df_res.set_index(['Cycle','ID'])
    df_res[['kout','kin','b','c','delta_kout','delta_kin','delta_b','delta_c','ChiSq']]=df_res[['kout','kin','b','c','delta_kout','delta_kin','delta_b','delta_c','ChiSq']].apply(pd.to_numeric)
    
    return df_res



print("Fitting model.....")
df_res = fit_LEXY_model(df_data)

# merge with prelit average
grp = df_data.groupby(['Cycle','ID'])
prelit_avg = grp[['meanint_LEXY','eccen','meanint_nucl','area','perimeter']].mean().reset_index().set_index(['Cycle','ID'])
df_res = df_res.join(prelit_avg)
df_res.head()

# save results
df_res.to_csv(analysis_dir+'model_fit_results.csv')



# #check the results
# fig = plt.figure(figsize=(8,7))
# #for i in np.arange(df_res.shape[0]):
# for i in np.arange(3):
#     row_res = df_res.iloc[i]
#     name = row_res.name
#     popt = tuple(row_res[['kout','kin','b','c']].values)

#     df_select = df_data.loc[(df_data['Cycle']==name[0]) & (df_data['ID']==name[1])]
    
#     x = df_select['macro_T']
#     y = df_select['meanint_LEXY_normed']
#     dy = df_select['stdint_LEXY_normed']/np.sqrt(df_select['area'])  #sem
    
#     plt.scatter(x,y,alpha=0.6)

#     xx = np.linspace(x.min(),x.max(),num=100)
#     yy = LEXY_model(xx,*popt)
#     plt.plot(xx,yy,'k')

