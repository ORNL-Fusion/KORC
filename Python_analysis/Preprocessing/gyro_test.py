import numpy as np
import matplotlib.pyplot as plt
import h5py

qe=1.6022e-19
me=9.1094e-31
mp = 1.6726e-27 #% Proton mass (kg)
c=2.9979e8
mu0=4*np.pi*10**(-7)

#%%loading data from KORC output file
#filename=r"/home/21b/KORC_RUNS/TEST0/OUT_NVHPC_ACC/file_0.h5"
filename=r"/home/21b/KORC/test/solovev/rank_1/file_0.h5"

iilen=102

time=np.zeros(iilen)
Kyy=np.zeros(iilen)
Kzz=np.zeros(iilen)

with h5py.File(filename,'r',locking=False) as f:
    for ii in range(0,iilen):
        time[ii]=f[str(ii)]['time'][0]
        Kyy[ii]=f[str(ii)]['spp_1']['X'][1][0]
        Kzz[ii]=f[str(ii)]['spp_1']['X'][2][0]
        

#%%calculate analytical orbit

sign=1

B=1.

K=10e6
g=K*qe/(me*c**2)+1
v=c*np.sqrt(1-1/g**2)

omega=qe*B/(g*me)
xi=2*np.pi
eta=np.pi/2
etadeg=10
eta=np.deg2rad(etadeg)

rL=v*np.sin(eta)/omega

tgyro=2*np.pi/omega

x0=0
y0=0
z0=0

vx0=0
#vy0=-v*np.cos(np.pi/4)
#vz0=v*np.cos(np.pi/4)
vy0=v
vz0=0

tsteps=209
#tsteps=14
dt=4.8e-13
#dt=0.7142857143E-11

tt=np.linspace(0,dt*tsteps,tsteps+1)

Ayy=sign*(vy0/omega*np.sin(omega*tt)-vz0/omega*(1-np.cos(omega*tt)))
Azz=sign*(vz0/omega*np.sin(omega*tt)+vy0/omega*(1-np.cos(omega*tt)))

#%% Plotting

plotyzorbit=0
plotyorbit=1
        
if plotyzorbit==1:
    fig,ax=plt.subplots()
    ax.plot(Kyy,Kzz,color='b')
    ax.plot(Ayy,Azz,color='r')
    ax.plot(Ayy[0],Azz[0],'o',color='k')
    ax.plot(Ayy[137],Azz[137],'*',color='k')
    ax.plot(Kyy[0],Kzz[0],'o',color='k')
    
    ax.set(xlabel='y (m)', ylabel='z (m)',
           title='y-z orbit')
    
    ax.legend(['KORC','Analytic'])
    
    ax.grid()
    plt.gca().set_aspect('equal')
    
    plt.show()

if plotyorbit==1:
    fig,ax=plt.subplots(2)
    ax[0].plot(time,Kyy,color='b')
    ax[0].plot(tt,Ayy,color='r')
    ax[0].plot(tt,(Kyy-Ayy)/Ayy,color='k')

    ax[0].set_ylim([-.1,-1])

    ax[0].set(xlabel='t (s)', ylabel='y (m)',
           title='y orbit')
    
    ax[0].legend(['KORC','Analytic'])
    
    ax[0].grid()
    
    ax[1].plot(time,Kzz,color='b')
    ax[1].plot(tt,Azz,color='r')
    ax[1].plot(tt,(Kzz-Azz)/Azz,color='k')
    
    ax[1].set(xlabel='t (s)', ylabel='z (m)',
           title='y orbit')
    
    ax[1].legend(['KORC','Analytic'])
    
    ax[1].grid()
    
    plt.show()
