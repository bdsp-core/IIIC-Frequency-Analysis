from fcn_startsetup import *
from fooof import FOOOF
from joblib import Parallel, delayed
from mne.filter import notch_filter,filter_data
import hdf5storage as hs
import scipy.io as sio
import numpy as np
import pandas as pd
import numpy.matlib
import math
import os

# global var
bipolar_channels=['FP1-F7','F7-T3','T3-T5','T5-O1','FP2-F8','F8-T4','T4-T6','T6-O2','FP1-F3','F3-C3','C3-P3','P3-O1','FP2-F4','F4-C4','C4-P4','P4-O2','FZ-CZ','CZ-PZ']
mono_channels=['FP1','F3','C3','P3','F7','T3','T5','O1','FZ','CZ','PZ','FP2','F4','C4','P4','F8','T4','T6','O2','EKG']
fooofgap=np.empty((1,5))
fooofgap[:]=np.nan
fooofgap=fooofgap[0]
thr_bw=0.5
freq_range=[.5,3]

# callbacks 
def fcn_computeSpectra(x,Fs):
    N=len(x)
    xdft=np.fft.fft(x)
    xdft=xdft[0:int(N/2)]
    psdx=(1/(Fs*N))*np.square(abs(xdft))
    psdx[1:-2]=2*psdx[1:-2]
    return psdx

def fcn_getBanana(X):
    bipolar_ids = np.array([[mono_channels.index(bc.split('-')[0]),mono_channels.index(bc.split('-')[1])] for bc in bipolar_channels])
    bipolar_data = X[bipolar_ids[:,0]]-X[bipolar_ids[:,1]]
    return bipolar_data

def fcn_rdafooof(seg,freqs, Fs):
    fooof=[]
    for k in range(0,len(seg)):
        s=fcn_computeSpectra(seg[k],Fs)
        spectrum=s[0]
        fm = FOOOF()
        freqs = np.array(freqs)
        spectrum = np.array(spectrum)
        try:
            fm.fit(freqs, spectrum, freq_range)
        except Exception as e:
            print(f"Error fitting FOOOF model for segment {k+1}: {e}")
            return 69420, 69420

        tmp=fm.peak_params_
        if len(tmp)==0:
            x=fooofgap
            fooof.append(x)
        else:
            x=np.append(np.append(tmp[0],fm.error_),fm.r_squared_)
            fooof.append(x)

    fooof=np.array(fooof)
    bw=np.round(100*fooof[:,2])/100 
    idx=np.where(bw<=thr_bw)
    fooof_=fooof[idx]
    channels_=np.array(bipolar_channels)[idx[0]]
    return fooof_,channels_

def process_file(sourceDir, scoreDir, targetDir, fname):
    res={"model_output":[],"model_prediction":[],"event_frequency":[],"bandwidth":[],"power":[],"fit_error":[],"r_squared":[],"spatial_extent":[],"channels":[]}
    csvpath = targetDir+fname.replace(".mat","_rdafreq.csv") 
    if os.path.isfile(csvpath):
        print('--alr done ' + fname)
        return
    
    # read raw EEG
    try:
        mat=hs.loadmat(sourceDir+fname)
    except Exception as ee:
        mat=sio.loadmat(sourceDir+fname)
        
    data=mat["data"]
    data=np.where(np.isnan(data),0,data)
    
    Fs=mat["Fs"]
    nWin=math.ceil(len(data[0])/(2*Fs))
    
    # filters to denoise
    data=notch_filter(data,Fs,60,n_jobs=-1,verbose="ERROR")
    data=filter_data(data,Fs,0.5,40,n_jobs=-1,verbose="ERROR")

    # L-bipolar
    data=fcn_getBanana(data)
    data=np.array(data)

    # read model output
    scrFile=scoreDir+"Data"+fname.replace(".mat","_score.csv")
    P_model=np.array(pd.read_csv(scrFile))
    P_model=np.concatenate((np.matlib.repmat(P_model[0],2,1),P_model),axis=0)
    if len(P_model)<nWin:
        P_model=np.concatenate((P_model,np.matlib.repmat(P_model[-1],nWin-len(P_model),1)),axis=0)  
    
    winsize=14*Fs
    stepsize=2*Fs 
    for i in range(0,nWin-7):
        eeg_start=int(i*stepsize)
        eeg_end=int(eeg_start+winsize) 
        
        seg=data[:16,eeg_start:eeg_end]
            
        scr_center=i+3
        yp=P_model[scr_center]
    
        res["model_output"].append(yp)
        res["model_prediction"].append(np.argmax(yp))
        
        if np.argmax(yp) in [4,5]: # RDAs
            N=len(seg[0])     
            freqs=np.arange(0,Fs/2,Fs/N) 
            fooof,channels=fcn_rdafooof(seg,freqs, Fs)
            print(fooof)
            print(channels)
            if np.any(fooof == 69420) or np.any(channels == 69420):
                res["event_frequency"].append(-2)
                res["power"].append(-2)
                res["bandwidth"].append(-2)
                res["fit_error"].append(-2)
                res["r_squared"].append(-2)
                res["spatial_extent"].append(-2)
                res["channels"].append(-2)
            else:
                if np.argmax(yp)==4: # LRDA
                    if len(channels)!=0:
                        scr=np.square(fooof[:,1])*fooof[:,4]/(fooof[:,2]*fooof[:,3]) 
                        idx=np.argsort(scr)                  
                        if len(idx)>1:
                            fooof=fooof[idx[[-1,-2]]]
                            channels=channels[idx[[-1,-2]]]
                                
                if len(channels)==0:
                    res["event_frequency"].append(-2)
                    res["power"].append(-2)
                    res["bandwidth"].append(-2)
                    res["fit_error"].append(-2)
                    res["r_squared"].append(-2)
                    res["spatial_extent"].append(-2)
                    res["channels"].append(-2)                 
                else:
                    res["event_frequency"].append(np.median(fooof[:,0]))
                    res["power"].append(np.median(fooof[:,1]))
                    res["bandwidth"].append(np.median(fooof[:,2]))
                    res["fit_error"].append(np.median(fooof[:,3]))
                    res["r_squared"].append(np.median(fooof[:,4]))
                    res["spatial_extent"].append(len(channels)/16)
                    res["channels"].append(channels)

        else: # Skip if not RDAs
            res["event_frequency"].append(-1)
            res["power"].append(-1)
            res["bandwidth"].append(-1)
            res["fit_error"].append(-1)
            res["r_squared"].append(-1)
            res["spatial_extent"].append(-1)
            res["channels"].append(-1)    

    res=pd.DataFrame(res)
    res.to_csv(targetDir+fname.replace(".mat","_rdafreq.csv"))
# main  
if __name__=='__main__':
    sourceDir="/bdsp/staging/TEEGLLTEEG/Datasets/continuousEEG/"
    scoreDir ="/pdDetector/finalfiles/"
    targetDir="/testing_pd/rdaFreqDetector/output2/"
    
    if not os.path.exists(targetDir):
        os.makedirs(targetDir)
        
    fnames=os.listdir(sourceDir)
    
    Parallel(n_jobs=32, verbose=5)(
        delayed(process_file)(sourceDir, scoreDir, targetDir, fname) for fname in fnames
    )
            
        