import numpy as np
import matplotlib.pyplot as plt

c = 2.99792458E8 # Speed of light (m/s)
qe = 1.60217662E-19 # Electron charge (C)
me = 9.10938356E-31 # Electron mass (kg)
kB = 1.3807E-23 # Boltzmann constant (J/K)
mu0 = (4E-7)*np.pi # Magnetic permeability [C/(m^2*s)]
ep0 = 8.854E-12 # Electric permitivity [C^2/(N*m^2)]
re = qe**2/(4*np.pi*ep0*me*c**2) # classical electron radius (m)
alpha_f = 1/137 # fine structure 

#%% calculating

B=2.2 #(T)
E=1 #(V/m)
EdotB=-1;

rad_damp_time=6*np.pi*ep0*(me*c)**3/(qe**4*B**2) #(s)

Klog=np.linspace(4,8,100) #(eV)
K=10**Klog #(eV)
g=1+K/(me*c**2/qe)
p=me*c*np.sqrt(g**2-1) #(kg m/s)

eta_deg=10
eta_rad=np.deg2rad(eta_deg)

dp_E=-EdotB*qe*E*np.cos(eta_rad)*np.ones(np.shape(K))
#dp_ECH=-EdotB*qe*E_CH*np.cos(eta_rad)
    
dp_R=-g*p*(1-np.cos(eta_rad)**2)/rad_damp_time
#dxi_R(:,ii)=cos(eta_test(ii)).*(1-cos(eta_test(ii)).^2)./(gamma_test*rad_damp_time)
#dxi_E(:,ii)=-JdotB*qe*E_test*(1-cos(eta_test(ii)).^2)./pmag_test

#%% plotting

K_vs_dp=1

if K_vs_dp==1:
    fig,ax=plt.subplots()
    
    
    ax.plot(K,dp_E,'r')
    ax.plot(K,-dp_R,'b')
    
    ax.set(xlabel='$K\,(\\mathrm{eV})$', ylabel='$dp/dt\,(\\mathrm{kg\cdot m/s^2})$')
    ax.grid()
    
    ax.legend(['dp_E','-dp_R'])
    
    plt.xscale('log')
    plt.yscale('log')
    
    #plt.axis([1.29e-2,1.3e-2,-400,400])

    ptitle=f"B={B}T,E={E}V/m,eta={eta_deg}deg"
    plt.title(ptitle)
    
    plt.savefig("K_vs_dp.pdf", format="pdf", bbox_inches="tight")
    plt.show()