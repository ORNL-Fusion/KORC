import numpy as np
import matplotlib.pyplot as plt
import h5py

qe=1.6022e-19
me=9.1094e-31
c=2.9979e8
mu0=4*np.pi*10**(-7)

#%%

K0=1E7 #in eV
K1=1E6 #in eV
K2=1E5 #in eV

g0=1+K0*qe/(me*c**2)
g1=1+K1*qe/(me*c**2)
g2=1+K2*qe/(me*c**2)

dt0=1e-2
dt1=dt0/g1*g0
dt2=dt0/g2*g0