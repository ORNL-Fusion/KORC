import numpy as np
import matplotlib.pyplot as plt
import h5py

qe=1.6022e-19
me=9.1094e-31
c=2.9979e8
mu0=4*np.pi*10**(-7)

#%%loading data from KORC output file

filename=r"/home/21b/KORC_RUNS/LOCAL/TEST1/OUT_GNU/simulation_parameters.h5"

with h5py.File(filename,'r') as f:
    num_snapshots=f['simulation']['num_snapshots'][0]+1  
    output_cadence=f['simulation']['output_cadence'][0]
    ppp=f['species']['ppp'][0]
       
filename=r"/home/21b/KORC_RUNS/LOCAL/TEST1/OUT_GNU/file_0.h5"

t_steps=np.zeros(num_snapshots,dtype=np.uint32) 
time=np.zeros((num_snapshots,ppp))
Kxx=np.zeros((num_snapshots,ppp))
Kyy=np.zeros((num_snapshots,ppp))
Kzz=np.zeros((num_snapshots,ppp))
Kvx=np.zeros((num_snapshots,ppp))
Kvy=np.zeros((num_snapshots,ppp))
Kvz=np.zeros((num_snapshots,ppp))
Kbx=np.zeros((num_snapshots,ppp))
Kby=np.zeros((num_snapshots,ppp))
Kbz=np.zeros((num_snapshots,ppp))
flagCon=np.zeros((num_snapshots,ppp))

with h5py.File(filename,'r') as f:
    for ii in range(0,num_snapshots):
        t_steps[ii]=output_cadence*ii
        time[ii]=f[str(t_steps[ii])]['time'][0]
        Kxx[ii]=f[str(t_steps[ii])]['spp_1']['X'][0][:]
        Kyy[ii]=f[str(t_steps[ii])]['spp_1']['X'][1][:]
        Kzz[ii]=f[str(t_steps[ii])]['spp_1']['X'][2][:]
        Kvx[ii]=f[str(t_steps[ii])]['spp_1']['V'][0][:]
        Kvy[ii]=f[str(t_steps[ii])]['spp_1']['V'][1][:]
        Kvz[ii]=f[str(t_steps[ii])]['spp_1']['V'][2][:]
        Kbx[ii]=f[str(t_steps[ii])]['spp_1']['B'][0][:]
        Kby[ii]=f[str(t_steps[ii])]['spp_1']['B'][1][:]
        Kbz[ii]=f[str(t_steps[ii])]['spp_1']['B'][2][:]
        flagCon[ii]=f[str(t_steps[ii])]['spp_1']['flagCon'][:]



#%% Plotting

plotyorbit=1       

if plotyorbit==1:
    fig,ax=plt.subplots(2)
    ax[0].plot(time,Kyy,color='b')


    ax[0].set(xlabel='t (s)', ylabel='y (m)',
           title='y orbit')
    
    ax[0].legend(['KORC','Analytic'])
    
    ax[0].grid()
    
    ax[1].plot(time,Kzz,color='b')

    
    ax[1].set(xlabel='t (s)', ylabel='z (m)',
           title='y orbit')
    
    ax[1].grid()
    
    plt.show()
