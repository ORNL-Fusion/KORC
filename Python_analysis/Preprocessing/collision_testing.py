import numpy as np
from scipy.interpolate import RegularGridInterpolator
import scipy as sci
import os
import matplotlib.pyplot as plt

#%% input variable
c_imp=0
E_c=5.0
nAr1_c = 0
nD1_c = 100

#% parameter setting
NPinit=100
tmax=50*1e-03  #%[s]
#% we consider mono-velocity pitch and radius as initial condition of REs
rho_ini = [0,1]   #% normalized radius (from 0 to 1, where rho = 1 corresponds to radius to minor radius Ra)
KEmag_ini = 5*1e6 #%(eV)  default is mono 2 Mev
Xi_ini = np.cos(0.174533) #% (rad) default is mono 10 deg

#%%-----------3D SDE code-----------
#% read & save files
prob_filename = 'LA_prob0.npy.npz'
save_filename = 'Result1.npy.npz'

if os.path.isfile(prob_filename):   #% if we need to construct prob.mat. Prob file depends on min/max gamma
    calcmoller=0
else:
    calcmoller=1

#% we consider 4 impurity species, where first initial try fixed nAr0 and
#% nD0. Here we do not have any impurity profile function, we just assume it
#% is constant for the whole simulaton.
nAr1_c=nAr1_c * 1e19   #%[/m**3]
nD1_c=nD1_c * 1e19     #%[/m**3]

nAr0_c=0#1e19 #%[/m**3]
nD0_c=0#1e19 #%[/m**3]



#%CONSTANTS
c = 2.99792458E8 #% Speed of light (m/s)
qe = 1.60217662E-19 #% Electron charge (C)
me = 9.10938356E-31 #% Electron mass (kg)
ep0 = 8.854E-12 #% Electric permitivity [C**2/(N*m**2)]
re = qe**2/(4*np.pi*ep0*me*c**2) #% classical electron radius (m)


#%% plasma setting default
#%MODELS
Eaccel=True
#% logical flag for including electric field acceleration
Rad=True
#% logical flag for including radiation acceleration
pitch_angle=True
#% include physics of collisional pitch angle scattering
slowing_down=True
#% include physics of collisional slowing down
par_diffusion=False
#% include physics of collisional parallel momentum diffusion
conservative=True
#% choice of having Boltzmann large angle collision operator conserve
constantClog=True

#% B0 in physical units (T)
B0 = 39#3.388 #% [T]
#%FIELD
B=B0   #% in physical units (T)
R0=2.616 #% major radius
rad_damp_time=6*np.pi*ep0*(me*c)**3/(qe**4*B**2)
#% radiation damping time scale
JdotB=-1
#%relative direction of current and magnetic field
k=5
KE_min = 1e-02  #% M(ev)
KE_max = 20    #% Mev

gammin=1+KE_min*1e6*abs(qe)/(me*c**2)          
#gammin=1.079266836027482#(gammin+1)/2

pmagmin=me*c*np.sqrt(gammin**2-1)

gam_ini=1+KEmag_ini/(me*c**2)*qe
omega_c=qe*B/(gam_ini*me)
dt_re=2*np.pi/omega_c

#%% Initialization
#%PARTICLES and Particle initialization
NP=100*NPinit
#% number of REs
NP_subset = round(0.7*NP)

rho=np.array(np.zeros(NP))
Vmag=np.array(np.zeros(NP))
Xi=np.array(np.zeros(NP))

#% Electrons position initialziation
EE0mag_ini = KEmag_ini+(me*c**2)/abs(qe)
rho[0:NPinit]=np.random.rand(NPinit) * (rho_ini[1]-rho_ini[0])+rho_ini[0]
Vmag[0:NPinit]=c*np.sqrt(1-(1./(EE0mag_ini*abs(qe)/(me*c**2)))**2)
Xi[0:NPinit]=Xi_ini

#%
flagdiss=np.array(np.zeros(NP))
flagdiss[0:NPinit]=1
#%initializing dissipation RE array
flagSec=np.array(np.zeros(NP))
flagSec[0:NPinit]=1
#%initializing able secondary RE array
flagRE=np.array(np.zeros(NP))
flagRE[0:NPinit]=1


#% snapshot for save file
n_snapshot = 201
IRE_time = np.array(np.zeros(n_snapshot))
IRE_Vmag = np.array(np.zeros((n_snapshot,NP)))
IRE_Xi = np.array(np.zeros((n_snapshot,NP)))
IRE_flagdiss = np.array(np.zeros((n_snapshot,NP)))
IRE_NREevo = np.array(np.zeros(n_snapshot))
IRE_NREevoact = np.array(np.zeros(n_snapshot))

#% Impurity
#% normalized effective length scales for Ne, Ar, D
aNe=[111,100,90,80,71,62,52,40,24,23,0]
aAr=[96,90,84,78,72,65,59,53,47,44,41,38,25,32,27,21,13,13,0]
aH=[274,0]

#%mean excitation energy of ion
INe=[137.2,165.2,196.9,235.2,282.8,352.6,475.0,696.8,1409.2,1498.4,np.inf]
IAr=[188.5,219.4,253.8,293.4,339.1,394.5,463.4,568.0,728.0,795.9,879.8,
     989.9,1138.1,1369.5,1791.2,2497.0,4677.2,4838.2,np.inf]
IH=[14.9916,np.inf]

#% Input parameters for Optimization
num_impurity_species = 4   #% multiple impurity species
Zo = np.array(np.zeros(num_impurity_species))
Zj = np.array(np.zeros(num_impurity_species))
Zj[0] = 0 #% Ar0
Zo[0] = 18
Zj[1] = 1  #% Ar+1
Zo[1] = 18
Zj[2] = 0  #% H0 (D)
Zo[2] = 1
Zj[3] = 1  #% H+1
Zo[3] = 1

IZj = np.array(np.zeros(num_impurity_species))
aZj = np.array(np.zeros(num_impurity_species))
IZj[0] = IAr[0]*abs(qe)    #% Ar
aZj[0] = aAr[0]
IZj[1] = IAr[1]*abs(qe)   #% Ar+1
aZj[1] = aAr[1]
IZj[2] = IH[0]*abs(qe)     #% H (D)
aZj[2] = aH[0]
IZj[3] = IH[1]*abs(qe)     #% H+1
aZj[3] = aH[1]


#%% Before pushing, computing dt and collisions parameters:
#% Impurity for dt and LA prob calculation
E = E_c
Te = 2000                 #%(eV)
Te= Te*abs(qe)         #% temperature of background thermal electrons (J)

ntotal_0 = Zo[0]*nAr0_c + Zo[1]*nAr1_c + Zo[2]*nD0_c + Zo[3]*nD1_c
ne_0 = Zj[0]*nAr0_c + Zj[1]*nAr1_c + Zj[2]*nD0_c + Zj[3]*nD1_c
Zeff_0=(nAr1_c+nD1_c)/ne_0

nION_frac_0 = np.array(np.zeros((num_impurity_species,1)))
nION_frac_0[0] = nAr0_c / ne_0
nION_frac_0[1] = nAr1_c / ne_0
nION_frac_0[2] = nD0_c / ne_0
nION_frac_0[3] = nD1_c / ne_0


#% Coulomb logarithm only for determining dt and prob
vth=np.sqrt(2*Te/me)     #% electron thermal speed
#% electron thermal speed
gammavth=1/np.sqrt(1-(vth/c)**2)

if constantClog:
    Zeff_0=2
    Clog0_0=20
else:
    Clog0_0=14.9-np.log(ne_0/10**(20))/2+np.log(Te/(10**(3)*abs(qe))) #%Coloumb logarithm
gammac_0=ne_0*abs(qe**4)*Clog0_0/(4*np.pi*ep0**2) #%collisional gamma
tau_c0=me**2*c**3/gammac_0
rel_col_freq0=1/tau_c0
#%collision frequency based off of relativistic velocity
E_CH0=rel_col_freq0*me*c/abs(qe)
if constantClog:
    Clog_ee0=Clog0_0
    Clog_ei0=Clog0_0
else:
    Clog_ee0=Clog0_0+np.log(1+(2*(gammin-1)*(c/vth)**2)**(k/2))/k
    Clog_ei0=Clog0_0+np.log(1+(2*pmagmin/(me*vth))**k)/k
E_CH0=E_CH0*(Clog_ee0/Clog0_0)*(ntotal_0/ne_0)

p_crit=me*c/np.sqrt(E/E_CH0-1)
gam_crit=np.sqrt((p_crit/(me*c))**2+1)
KE_crit=(gam_crit-1)*me*c**2/qe

#%% LA initial only dependes on gammin,gammax
neta0=100
neta1=101
ngam0=102
ngam1=103

#% setting variable arrays
eta0=np.linspace(0,np.pi,neta0)
eta1=np.linspace(0,np.pi,neta1)

gamoffset=gammin
gammax=20.569511968719837

gam0=10**(10**(np.linspace(np.log10(np.log10(gamoffset)),np.log10(np.log10(gammax)),ngam0)))
gam1=10**(10**(np.linspace(np.log10(np.log10(gamoffset)),np.log10(np.log10((gammax+1)/2)),ngam1)))

KE0=(gam0-1)*me*c**2/qe
KE1=(gam1-1)*me*c**2/qe

pmag1=me*c*np.sqrt(gam1**2-1)

[gam00,gam11,eta00,eta11]=np.meshgrid(gam0,gam1,eta0,eta1,indexing='ij')
[gam000,eta000]=np.meshgrid(gam0,eta0,indexing='ij')

#% prob depends on:
if calcmoller == 1:
    p00=me*c*np.sqrt(gam00**2-1)
    p11=me*c*np.sqrt(gam11**2-1)

    #%pitch deflection
    cosgam=np.sqrt((gam00+1)*(gam11-1)/((gam00-1)*(gam11+1)))

    sinsq=np.sin(eta00)**2*np.sin(eta11)**2
    cossq=np.cos(eta00)*np.cos(eta11)
    cossq1=(cosgam-cossq)**2

    pitchprob=1/(np.pi*np.sqrt(sinsq-cossq1))

    pitchprob1=pitchprob
    pitchprob1[np.imag(pitchprob1)!=0]=0

    #%Moller scattering differential cross-section

    dsigdgam=2*np.pi*re**2/(gam00**2-1)*((gam00-1)**2*gam00**2/
        ((gam11-1)**2*(gam00-gam11)**2)-(2*gam00**2+2*gam00-1)/
        ((gam11-1)*(gam00-gam11))+1)

    #%Boltzmann collision operator (Eq1 in Matt's paper)

    S_CB=1/(2*np.pi*me**3*c**2)/(p11*gam11)*(p00/gam00)*pitchprob1*dsigdgam
    S_CB[np.isnan(S_CB)]=0.

    secthreshgam=np.ones(np.shape(gam11))
    secthreshgam[gam11>(gam00+1)/2]=0
    #%secthreshgam(gam00<gamcrit)=0

    S_CB=S_CB*secthreshgam
    
    p2S_CB=p11**2*S_CB

    #% Number of secondary REs generated by an incident RE per unit time
    intp1=np.trapz(p2S_CB,x=pmag1,axis=1)

    eta110=eta11[:,1,:,:]

    prob=2*np.pi*np.trapz(np.sin(eta110)*intp1,x=eta1,axis=2)
    #prob[np.isnan(prob)]=0.

    #% maximum CB of secondary RE as a function of primary RE
    maxCB=np.array(np.zeros((ngam0,neta0)))

    for ii in range(0,neta0):
        for jj in range(0,ngam0):
            temp=p2S_CB[jj,:,ii,:]

            maxCB[jj,ii]=np.max(temp)
        
    
    np.savez(prob_filename,mat1=maxCB,mat2=prob)
else:
    loaded_data=np.load(prob_filename)
    maxCB=loaded_data['mat1']
    prob=loaded_data['mat2']


maxCB = maxCB*ntotal_0
prob = prob*ntotal_0

#%% plotting for LAC

plotdata=1
plotcdf=1

if plotcdf==1:
    
    if neta0==100:
        gam0ind=81
        eta0ind=5
    elif neta0==50:
        gam0ind=51
        eta0ind=3
    
    [pm11,sineta11]=np.meshgrid(pmag1,np.sin(eta1),indexing='ij')
    
    dpmag1=np.diff(pmag1,append=0)
    dpmag1[-1]=dpmag1[-2]
    
    [dpm11,sineta11]=np.meshgrid(dpmag1,np.sin(eta1),indexing='ij')
    
    tmp=dpm11*pm11**2*sineta11*S_CB[gam0ind,:,eta0ind,:]
    tmp1=tmp
    tmp1=tmp1.flatten()
    
    tmp2=np.cumsum(tmp1)
    tmp2=tmp2/np.max(tmp2)

    xrange=np.linspace(0,1,tmp2.shape[0])
    
    fig,ax=plt.subplots()
    
    ax.plot(xrange,tmp2,linewidth=2,color='k')
    
    ax.plot([0.27,1.05],[0.753,0.753],linestyle='--',color='k')
    ax.plot([0.27,0.27],[-.05,0.753],linestyle='--',color='k')
    
    ax.set(xlabel='$\mathcal{X}(p,\eta)$', ylabel='CDF')
    ax.set(title='$\\mathcal{K}_0=10\,{\\mathrm{MeV}}$, $\eta_0=11^\circ$')
    ax.grid()
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(-0.05, 1.05)
    
    ax2 = ax.twinx()
    ax2.set(ylabel='$\\mathcal{U}(0,1)$')
    ax2.set_ylim(-0.05, 1.05)
    
    #ax.set_xscale('log')
    
    
    plt.savefig("CDF.png", format="png", bbox_inches="tight")
    plt.show()  


if plotdata==1:
    fig,ax=plt.subplots()
    
    if neta0==100:
    
        gam0ind=86
        eta0ind=5
        
        gam1ind=86
        eta1ind=5
        
    elif neta0==50:
    
        gam0ind=50
        eta0ind=3
        
        gam1ind=44
        eta1ind=3
    
    #tmp=prob
    #tmp=maxCB
    #tmp=p11[gam0ind,:,eta0ind,:]
    #tmp=pitchprob[gam0ind,:,eta0ind,:]
    #tmp=pitchprob1[gam0ind,:,eta0ind,:]
    #tmp=dsigdgam[gam0ind,:,eta0ind,:]
    tmp=S_CB[gam0ind,:,eta0ind,:]
    #tmp=p2S_CB[gam0ind,:,eta0ind,:]
    #tmp=secthreshgam[gam0ind,:,eta0ind,:]
    #tmp=p11**2*S_CB
    #tmp=tmp[gam0ind,:,eta0ind,:]
    #tmp=intp1[:,:,eta1ind]
    #tmp=np.sin(eta110[:,:,eta1ind])
    #tmp=np.sin(eta110)*intp1
    #tmp=tmp[:,:,eta1ind]
    #tmp=eta110[:,:,eta1ind]
    
    # Create a modified colormap with 'white' as the "under" color
    cmap = plt.cm.viridis
    cmap_modified = cmap.with_extremes(under='white')
    
    plotlog=1
    
    if plotlog==1:
        #ct=ax.pcolormesh(KE1,eta1,np.log10(np.transpose(tmp)),cmap=cmap_modified,vmin=.001)
        ct=ax.pcolormesh(KE1,eta1,np.log10(np.transpose(tmp)),cmap=cmap_modified,vmin=.001)
        
        cb=plt.colorbar(ct)
        cb.set_label('${\\rm log}_{10}(f_{RE})$', rotation=90)
    else:
        ct=ax.pcolormesh(KE1,eta1,np.transpose(tmp),cmap=cmap_modified,vmin=.001)
        #ct=ax.pcolormesh(KE0,eta0,np.transpose(tmp),cmap=cmap_modified,vmin=.001)
        
        cb=plt.colorbar(ct)
        cb.set_label('$f_{RE}$', rotation=90)
    
    
    
    ax.set(xlabel='$\\mathcal{K} (\\mathrm{eV})$', ylabel='$\eta (\\mathrm{rad})$')
    ax.grid()
    
    ax.set_xscale('log')
    
    
    plt.savefig("LAC.pdf", format="pdf", bbox_inches="tight")
    plt.show()  


#%% TIMINGS
KEmin=me*c**2*(gammin-1)/qe
EEmin=KEmin+(me*c**2)/np.abs(qe)
Vmagmin=c*np.sqrt(1-(1/(EEmin*np.abs(qe)/(me*c**2)))**2)
x=Vmagmin/vth

if constantClog:
    Clog_ee=Clog0_0
    Clog_ei=Clog0_0
else:
    Clog_ee=Clog0_0+np.log(1+(2*(gammin-1)*(c/vth)**2)**(k/2))/k
    Clog_ei=Clog0_0+np.log(1+(2*pmagmin/(me*vth))**k)/k

Chand=(sci.special.erf(x)-x*(2/np.sqrt(np.pi))*np.exp(-(x**2)))/(2*x**2)

CF = gammac_0*(Clog_ee/Clog0_0)*Chand/Te
CF_temp = CF
for ii in range(0,num_impurity_species):
    h_PS=(pmagmin/(me*c))*np.sqrt(gammin-1)*(me*c**2/IZj[ii])
    CF_temp = CF_temp + CF * nION_frac_0[ii] * (Zo[ii] - Zj[ii])/Clog_ee*(np.log(1+h_PS**k)/k-(Vmagmin/c)**2)

CF = CF_temp

CB_ei = gammac_0*(Clog_ee/Clog0_0)/(2*Vmagmin)*(Clog_ei/Clog_ee)*Zeff_0
CB_ei_temp = CB_ei
for ii in range(0,num_impurity_species):

    g_PS=2/3*(Zo[ii]**2-Zj[ii]**2)*np.log(((pmagmin/(me*c))*aZj[ii])**(3/2)+1)-2/3*(Zo[ii]-Zj[ii])**2*((pmagmin/(me*c))*aZj[ii])**(3/2)/(((pmagmin/(me*c))*aZj[ii])**(3/2)+1)

    CB_ei_temp=CB_ei_temp + CB_ei*nION_frac_0[ii]/Zeff_0/Clog_ei*g_PS

CB_ei = CB_ei_temp
CB_ee = gammac_0*(Clog_ee/Clog0_0) /(2*Vmagmin) * (sci.special.erf(x)-Chand+(vth*Vmagmin/c**2)**2/2)
CB = CB_ei + CB_ee


CA=gammac_0*(Clog_ee/Clog0_0)*Chand/(Vmagmin)
dCA=gammac_0*(Clog_ee/Clog0_0)/(gammin**3*me*Vmagmin**2)*((2*(gammin*Vmagmin/c)**2-1)*Chand+x*(2/np.sqrt(np.pi))*np.exp(-(x**2)))

nuS=CF/pmagmin
nuD=2*CB/pmagmin**2
nuP=2*CA/pmagmin**2
nudP=dCA/pmagmin

tau=min(1/nuS[0],1/nuD[0])
tau=min(tau,1/nuP)

#% TIMING steps
tau=min(tau,1/np.max(prob))
dt=tau/np.round(20)
numits=int(np.round(tmax/dt,-1))


Etime=np.array(np.zeros(numits+1))
NREevoact=np.array(np.zeros(numits+1))
NREevo=np.array(np.zeros(numits+1))

inds=np.where(flagRE>0)
lastind=inds[0][-1]

Etime[0]=0.
NREevoact[0]=lastind+1
NREevo[0]=lastind+1

pp=0

#%% Iterator
for ii in range(0,numits):
    #% R0 = profilesM_t.R0
    
    lastind=np.where(flagRE==1)
    lastind=lastind[0][-1]
    
    NRE=lastind+1
    N_imp = 1+np.tanh(c_imp*(rho[0:NRE]-0.5))
    N_imp=np.array(N_imp)

    nAr0 = nAr0_c*N_imp
    nAr1 = nAr1_c*N_imp
    nD0 = nD0_c*N_imp
    nD1 = nD1_c*N_imp

    ne=nAr1+nD1
    Zeff=(nAr1+nD1)/ne

    nION_frac = np.array(np.zeros((num_impurity_species,np.shape(ne)[0])))
    nION_frac[0,:] = nAr0/ne
    nION_frac[1,:] = nAr1/ne
    nION_frac[2,:] = nD0/ne
    nION_frac[3,:] = nD1/ne

    n_total = ne
    for ii_imp in range(0,num_impurity_species):
        n_total = n_total + ne*(Zo[ii_imp]-Zj[ii_imp])*nION_frac[ii_imp,:]

    #% relativistic factor of thermal electron
    if constantClog:
        Clog0=20
        Zeff=2
    else:
        Clog0=14.9-np.log(ne/10**(20))/2+np.log(Te/(10**(3)*abs(qe)))
    #%Coloumb logarithm
    gammac=ne*abs(qe**4)*Clog0/(4*np.pi*ep0**2)
    #%collisional gamma
    th_col_freq=gammac/(me**2*vth**3)
    #%collision frequency based off of electron thermal velocity
    rel_col_freq=gammac/(me**2*c**3)
    #%collision frequency based off of relativistic velocity
    E_CH=rel_col_freq*me*c/abs(qe)
    if constantClog:
        Clog_ee0=Clog0
        Clog_ei0=Clog0
    else:
        Clog_ee0=Clog0+np.log(1+(2*(gammin-1)*(c/vth)**2)**(k/2))/k
        Clog_ei0=Clog0+np.log(1+(2*pmagmin/(me*vth))**k)/k
    E_CH=E_CH*(Clog_ee0/Clog0)*(n_total/ne)
    #%Connor-Hastie electric field
    
    if np.mod(ii,numits/10)==0:
        X=['Iteration: ',str(ii)]
        print(X)

    gamma=1/np.sqrt(1-(Vmag[0:NRE]/c)**2)
    pmag=me*gamma*Vmag[0:NRE]

    #% Deterministic acceleration
    dpdet_E=np.array(np.zeros(NRE))
    dxidet_E=np.array(np.zeros(NRE))
    dpdet_R=np.array(np.zeros(NRE))
    dxidet_R=np.array(np.zeros(NRE))

    if Eaccel:
        dpdet_E=dpdet_E-JdotB*qe*E*Xi[0:NRE]*dt
        dxidet_E=dxidet_E-JdotB*qe*E*(1-Xi[0:NRE]**2)/pmag*dt

    if Rad:
        dpdet_R=dpdet_R-gamma*pmag*(1-Xi[0:NRE]**2)/rad_damp_time*dt
        dxidet_R=dxidet_R+Xi[0:NRE]*(1-Xi[0:NRE]**2)/(gamma*rad_damp_time)*dt

    dp_D=dpdet_E+dpdet_R
    dxi_D=dxidet_E+dxidet_R


    pmag=pmag+flagdiss[0:NRE]*dp_D
    Xi[0:NRE]=Xi[0:NRE]+flagdiss[0:NRE]*dxi_D

    gamma=np.sqrt(1+(pmag/(me*c))**2)
    Vmag[0:NRE]=pmag/(me*gamma)

    if sum(pmag<0) >0:
        print('pmag wrong')
        exit()

        
    Xi=[1-np.mod(Xii,1) if Xii>1 else Xii for Xii in Xi]
    Xi=[-1-np.mod(Xii,-1) if Xii<-1 else Xii for Xii in Xi]
    Xi=np.array(Xi)


    #% Stochastic Collision Operator
    dWV=np.sqrt(3*dt)*(-1+2*np.random.rand(NRE))
    dWT=np.sqrt(3*dt)*(-1+2*np.random.rand(NRE))
    x=Vmag[0:NRE]/vth

    if constantClog:
        Clog_ee=Clog0
        Clog_ei=Clog0
    else:
        Clog_ee=Clog0+np.log(1+(2*(gamma[0:NRE]-1)*(c/vth)**2)**(k/2))/k
        Clog_ei=Clog0+np.log(1+(2*pmag[0:NRE]/(me*vth))**k)/k

    Chand=(sci.special.erf(x)-x*(2/np.sqrt(np.pi))*np.exp(-(x**2)))/(2*x**2)

    CF = gammac*(Clog_ee/Clog0)*Chand/Te
    CF_temp = CF
    
    for jj in range(0,num_impurity_species):
        h_PS=(pmag/(me*c))*np.sqrt(gamma-1)*(me*c**2/IZj[jj])
        CF_temp = CF_temp + CF*nION_frac[jj,:]*(Zo[jj] - Zj[jj])/Clog_ee* \
            (np.log(1+h_PS**k)/k-(Vmag[0:NRE]/c)**2)

    CF = CF_temp

    CB_ei = gammac*(Clog_ee/Clog0)/(2*Vmag[0:NRE])*(Clog_ei/Clog_ee)*Zeff
    CB_ei_temp = CB_ei
    
    for jj in range(0,num_impurity_species):
        g_PS=2/3*(Zo[jj]**2-Zj[jj]**2)*np.log(((pmag/(me*c))*aZj[jj])**(3/2)+1)- \
            2/3*(Zo[jj]-Zj[jj])**2*((pmag/(me*c))*aZj[jj])**(3/2)/(((pmag/(me*c))*aZj[jj])**(3/2)+1)
        CB_ei_temp=CB_ei_temp + CB_ei*nION_frac[jj,:]/Zeff/Clog_ei*g_PS

    CB_ei = CB_ei_temp
    CB_ee = gammac*(Clog_ee/Clog0)/(2*Vmag[0:NRE])*(sci.special.erf(x)-Chand+(vth*Vmag[0:NRE]/c**2)**2/2)
    CB = CB_ei + CB_ee

    CA=gammac*(Clog_ee/Clog0)*Chand/(Vmag[0:NRE])

    dCA=gammac*(Clog_ee/Clog0)/(gamma**3*me*(Vmag[0:NRE])**2)* \
        ((2*(gamma*Vmag[0:NRE]/c)**2-1)*Chand+x*(2/np.sqrt(np.pi))*np.exp(-(x**2)))


    dpdet_C=np.array(np.zeros(NRE))
    dprand_C=np.array(np.zeros(NRE))

    if (slowing_down):
        dpdet_C=dpdet_C-CF*dt

    if (par_diffusion):
        dpdet_C=dpdet_C+dCA*dt
        dprand_C=dprand_C+np.sqrt(2*CA)*dWV

    dxidet_C=np.array(np.zeros(NRE))
    dxirand_C=np.array(np.zeros(NRE))

    if (pitch_angle):
        dxidet_C=dxidet_C-2*Xi[0:NRE]*CB/pmag**2*dt
        dxirand_C=dxirand_C-np.sqrt(2*CB)/pmag*np.sqrt(1-Xi[0:NRE]**2)*dWT


    dp_C=dpdet_C+dprand_C
    dxi_C=dxidet_C+dxirand_C

    pmag=pmag+flagdiss[0:NRE]*dp_C
    Xi[0:NRE]=Xi[0:NRE]+flagdiss[0:NRE]*dxi_C

    Xi=[1-np.mod(Xii,1) if Xii>1 else Xii for Xii in Xi]
    Xi=[-1-np.mod(Xii,-1) if Xii<-1 else Xii for Xii in Xi]
    Xi=np.array(Xi)

    gamma=np.sqrt(1+(pmag/(me*c))**2)
    Vmag[0:NRE]=pmag/(me*gamma)


    if max(gamma) > gammax:
        print('gamma out (over) of range')
        exit()

    #% tuning flagdiss   
    for jj in range(0,NRE):
        if gamma[jj]<=gammin*1.001:
            flagdiss[jj]=0
        if E<E_CH[jj]:
            flagSec[jj]=0

    #% large angle collision operator
    NREactive=sum(flagRE*flagdiss)
    ind0=np.where(flagRE*flagdiss*flagSec>0)[0]
    if ind0.size>0:
        
        newRE=np.shape(ind0)[0]
        
        prob0=np.random.rand(newRE)
        
        interp_prob = RegularGridInterpolator((gam0,eta0), prob, method='linear')
        prob_q=interp_prob([np.column_stack((gamma[ind0],np.acos(Xi[ind0])))])
       
        prob1 = dt*prob_q*n_total[ind0]/ntotal_0
        
        if (any(np.isnan(prob1)[0])):
            Y=['Min gamma: ',np.num2str(min(gamma[ind0]))]
            print(Y)
            print('NaN probability from Boltzmann source!!!')
            exit()

        if any((prob1>1)[0]):
            exit('multiple secondary REs generated in a time step!!!')
        
        ind0=np.delete(ind0,np.where(prob1<prob0))
        
    if ind0.size>0:
        
        newRE=np.shape(ind0)[0]
        
        flagRE[NRE:NRE+newRE]=1
        flagdiss[NRE:NRE+newRE]=1
        flagSec[NRE:NRE+newRE]=1

        rhoNew = rho[ind0]

        #% 2d pseudo-inverse cdf
        n_sample = 10000
        gamsecmax=(gamma[ind0]+1)/2
        gamsecmax = np.tile(gamsecmax, (n_sample,1))
        psecmax=me*c*np.sqrt(gamsecmax**2-1)

        gamma_ind0 = np.tile(gamma[ind0], (n_sample,1))
        pmag_ind0 = np.tile(pmag[ind0], (n_sample,1))

        ptrial=pmagmin+(psecmax-pmagmin)*np.random.rand(n_sample, newRE)

        gamtrial=np.sqrt(1+(ptrial/(me*c))**2)
        cosgam10=np.sqrt((gamma_ind0+1)*(gamtrial-1)/((gamma_ind0-1)*(gamtrial+1)))

        xitrial=-1+2*np.random.rand(n_sample, newRE)

        Xi_ind0 = np.tile(Xi[ind0], (n_sample,1))

        sinsq=(1-Xi_ind0**2)*(1-xitrial**2)
        cossq=Xi_ind0*xitrial
        cossq1=(cosgam10-cossq)**2
        pitchprob10=1/(np.pi*np.sqrt(sinsq-cossq1))
        
        pitchprob10[np.imag(pitchprob10)!=0]=0

        dsigdgam10=2*np.pi*re**2/(gamma_ind0**2-1)* \
            ((gamma_ind0-1)**2*gamma_ind0**2/ \
            ((gamtrial-1)**2*(gamma_ind0-gamtrial)**2)- \
            (2*gamma_ind0**2+2*gamma_ind0-1)/ \
            ((gamtrial-1)*(gamma_ind0-gamtrial))+1)

        disttrial= ptrial**2*n_total[ind0]/(2*np.pi*me**3*c**2)/ \
            (ptrial*gamtrial)*(pmag_ind0/ \
            gamma_ind0)*pitchprob10*dsigdgam10

        interp_maxS = RegularGridInterpolator((gam0,eta0), maxCB,method='linear')
        maxS_q=interp_maxS([np.column_stack((gamma[ind0],np.acos(Xi[ind0])))])

        maxS = np.tile(maxS_q,(n_sample,1))

        ratio = disttrial/maxS
        ratio[np.isnan(ratio)]=0.
        ratio = np.cumsum(ratio,axis=0)
        ratio = ratio/np.tile(ratio[-1,:], (n_sample,1))

        rand_temp = np.tile(np.random.rand(newRE), (n_sample,1))

        flag = rand_temp < ratio
        maxflag = np.argmax(flag,axis=0)

        ptrial = ptrial[maxflag,np.arange(flag.shape[1])]
        gamtrial = gamtrial[maxflag,np.arange(flag.shape[1])]
        xitrial = xitrial[maxflag,np.arange(flag.shape[1])]

        pmagNew=ptrial
        gamNew=gamtrial
        xiNew=xitrial

        if conservative:
            gammap0=gamma[ind0]
            #% gamma(ind0)=gamma(ind0)-gamtrial+gammavth(ind0)
            gamma[ind0]=gamma[ind0]-gamtrial+gammavth
            Vmag[ind0]=c*np.sqrt(1-1/gamma[ind0]**2)
            Xi[ind0]=(np.sqrt(gammap0**2-1)*Xi[ind0]- \
                np.sqrt(gamtrial**2-1)*xitrial)/np.sqrt(gamma[ind0]**2-1)
            #%error('manipulate state')

        Vmag[NRE:NRE+newRE]=pmagNew/(me*gamNew)
        Xi[NRE:NRE+newRE]=xiNew
        rho[NRE:NRE+newRE]=rhoNew

        Etime[ii]=ii*dt
        NREevo[ii]=NRE
        NREevoact[ii]=NREactive

        reducenum=0
        if reducenum==1:
            if (ii<=np.round(0.95*numits) and NRE>=NP):
                subset_idx = np.random.permutation(NP)[0:NP_subset]
                Vmag[0:NP_subset] = Vmag[subset_idx]
                Vmag[NP_subset:-1] = 0
    
                Xi[0:NP_subset] = Xi[subset_idx]
                Xi[NP_subset:-1] = 0
    
                rho[0:NP_subset] = rho[subset_idx]
                rho[NP_subset:-1] = 0
    
                flagdiss[0:NP_subset] = flagdiss[subset_idx]
                flagdiss[NP_subset:-1] = 0
    
                flagSec[0:NP_subset] = flagSec[subset_idx]
                flagSec[NP_subset:-1] = 0
    
                flagRE[0:NP_subset] = flagRE[subset_idx]
                flagRE[NP_subset:-1] = 0
    
            elif (ii>np.round(0.95*numits) and NRE>=NP):
                break

    else:
        Etime[ii]=ii*dt
        NREevo[ii]=NRE
        NREevoact[ii]=NREactive


    if np.mod(ii,np.round(numits/(n_snapshot-1)))==0:
        IRE_time[pp] = Etime[ii]
        IRE_Vmag[pp,:] = Vmag[0:NP]
        IRE_Xi[pp,:] = Xi[0:NP]
        IRE_flagdiss[pp,:] = flagdiss[0:NP]
        IRE_NREevo[pp] = NREevo[ii]
        IRE_NREevoact[pp] = NREevoact[ii]
        pp=pp+1


#%% Saving output data

ifsave=1

if ifsave==1:
    np.savez(save_filename,mat1=Etime,mat2=NREevo,mat3=NREevoact,mat4=IRE_flagdiss,mat5=IRE_time,mat6=IRE_Vmag,mat7=IRE_Xi,mat8=IRE_NREevo,mat9=IRE_NREevoact)

#%% plotting

plt.rcParams['contour.negative_linestyle'] = 'solid'

plotevo=1

if plotevo==1:
    
    tmp_time=IRE_time.copy()
    tmp_NREtot=IRE_NREevo.copy()
    tmp_NREact=IRE_NREevoact.copy()
    
    tmp_time=np.trim_zeros(tmp_time,trim='b')
    tmp_NREtot=np.trim_zeros(tmp_NREtot,trim='b')
    tmp_NREact=np.trim_zeros(tmp_NREact,trim='b')
    
    iffit=1
    
    if iffit==1:
        indrange0=int(np.round(np.shape(tmp_time)[0]*.5))
        
        def exp_func(x,b,c):
            return b*np.exp(c*x)
        
        popt_tot, pcov = sci.optimize.curve_fit(exp_func, tmp_time[indrange0:-1], tmp_NREtot[indrange0:-1], p0=[tmp_NREtot[0],1])
        popt_act, pcov = sci.optimize.curve_fit(exp_func, tmp_time[indrange0:-1], tmp_NREact[indrange0:-1], p0=[tmp_NREact[0],1])
        
        b_fit_tot,c_fit_tot=popt_tot
        b_fit_act,c_fit_act=popt_act
    
    fig,ax=plt.subplots()
    
    ax.plot(tmp_time,tmp_NREtot,'o-')
    ax.plot(tmp_time,tmp_NREact,'o-')
    
    if iffit==1:
        ax.plot(tmp_time[indrange0:-1],b_fit_tot*np.exp(c_fit_tot*tmp_time[indrange0:-1]))
        ax.plot(tmp_time[indrange0:-1],b_fit_act*np.exp(c_fit_act*tmp_time[indrange0:-1]))
    
    #ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.text(0.05,.85,r"$E/E_{\rm CH}=$"+str(np.round(E/E_CH[0],2))+"\n"+r"$\gamma \tau_c=$"+str(np.round(c_fit_act*tau_c0,3)),transform=ax.transAxes)
    
    ax.set(xlabel='$t\,({\\rm s})$', ylabel='$N_{RE}$')
    ax.grid()
    
    plt.savefig("NREevo.png", format="png", bbox_inches="tight")
    plt.show()






