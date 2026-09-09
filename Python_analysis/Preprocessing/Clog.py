import numpy as np

#%%

ne=10**(16) #in m^-3
Te=1*10**2  #in eV

Cbrag_lt50=23.4-1.15*np.log10(ne/(10**6))+3.45*np.log10(Te) 
Cbrag_gt50=25.3-1.15*np.log10(ne/(10**6))+2.30*np.log10(Te) 

Cwes=17.3-0.5*np.log(ne/(10**20))+1.5*np.log(Te/(10**3))
