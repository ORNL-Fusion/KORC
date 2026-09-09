import numpy as np
import matplotlib.pyplot as plt
import scipy

mu0=4*np.pi*10**(-7)
c = 2.99792458E8 #% Speed of light (m/s)
qe = 1.60217662E-19 #% Electron charge (C)
me = 9.10938356E-31 #% Electron mass (kg)
ep0 = 8.854E-12 #% Electric permitivity [C**2/(N*m**2)]
re = qe**2/(4*np.pi*ep0*me*c**2) #% classical electron radius (m)

#%% fitting varying elongation

radii = np.array([
    0.06591262, 0.09884386, 0.13175743, 0.16459691, 0.19729905, 0.22977723,
    0.26193569, 0.29365476, 0.32480509, 0.3552533 , 0.38487346, 0.41351341,
    0.44098735, 0.46713003, 0.49172032, 0.51454623, 0.5349404
])


kappa_r = np.array([
    1.22291, 1.22609, 1.2297 , 1.23414, 1.23984, 1.24711, 1.25581, 1.26617,
    1.27842, 1.29304, 1.31025, 1.33048, 1.35444, 1.38288, 1.417  , 1.45834,
    1.51087
])

initial_guess=[1,1,1]

def variation(radius,kap0,alpha,beta):
    
    var=kap0*(1+(radius/alpha)**beta)
    return var

popt, pcov = scipy.optimize.curve_fit(variation,radii,kappa_r,p0=initial_guess)

#%%calculating domains and fields

elonggeo=0
if elonggeo==1:
    
    R0=1.7
    Z0=0
    a=0.55

    Rw=1.9
    Zh=1.9

    R=np.linspace(R0-Rw/2,R0+Rw/2,1000)
    Z=np.linspace(Z0-Zh/2,Z0+Zh/2,1000)

    RR, ZZ = np.meshgrid(R, Z)
    
    #psib=1
    #a=1
    #kappa=1.5
    #psi=psib/a**2*(((RR**2-R0**2)**2/(4*R0**2)+ZZ**2/kappa**2))
    
    lam=0.284677
    B0=1.6
    q0=1.01
    kappa=1.5
    tri=0
    
    qa=q0*(1+(a/lam)**2)
    rm=np.sqrt((RR-R0)**2+(ZZ-Z0)**2)
    
    rm_elong=np.sqrt((RR-R0)**2+(ZZ-Z0)**2/kappa**2)
    
    q=q0*(1+(rm_elong/lam)**2)
    
    theta=np.atan2((ZZ-Z0),(RR-R0))
    psi=lam**2*B0/(2*q0)*np.log(1+(rm_elong/lam)**2)
    
    limR=R0+a*np.cos(np.linspace(0,2*np.pi,100))
    limZ=Z0+a*np.sin(np.linspace(0,2*np.pi,100))*kappa
    
    
    A_turb    = 1.620386820342208e-06
    mu_turb   = 0.32309683251166604
    sigma_turb= 0.017691893772250868
    
    g_r = A_turb * np.exp(-0.5 * ((rm_elong - mu_turb)/sigma_turb)**2)
    ballooning = 0.25 * (1.0 + np.cos(theta))**2
    dBr_norm_squared=g_r*ballooning
    
    Br = np.sqrt(dBr_norm_squared)
    
    RRdiffusion=1
    if RRdiffusion==1:
        KE=1e7*qe
        eta=np.deg2rad(10)
        
        gam=1+KE/(me*c**2)
        E=gam*me*c**2
        v=c*np.sqrt(1-(me*c**2/E)**2)
        vpll=v*np.cos(eta)

        Lc=np.pi*R0*q
        DRR=vpll*Lc*dBr_norm_squared
        
        ddeltaBnorm=g_r*ballooning*(-(rm_elong - mu_turb)/sigma_turb**2)
        dLc=np.pi*R0*q0*2*rm_elong/lam**2
        dDRR=vpll*(Lc*ddeltaBnorm+dLc*dBr_norm_squared)
        
        mu_drift=dDRR+DRR/rm_elong
        sigma=np.sqrt(2*DRR)
        
        nu_drift=mu_drift/rm_elong
        nu_RR=2*DRR/rm_elong**2
        
        q_mu=q0*(1+(mu_turb/lam)**2)
        dBr_norm_squared_mu = A_turb
        
        Lc_mu=np.pi*R0*q_mu
        DRR_mu=vpll*Lc_mu*dBr_norm_squared_mu
        
        ddeltaBnorm_mu=0
        dLc_mu=np.pi*R0*q0*2*mu_turb/lam**2
        dDRR_mu=vpll*(Lc_mu*ddeltaBnorm_mu+dLc_mu*dBr_norm_squared_mu)
        
        mu_drift_mu=dDRR_mu+DRR_mu/mu_turb
        sigma_mu=np.sqrt(2*DRR_mu)
        
        nu_drift_mu=mu_drift_mu/mu_turb
        nu_RR_mu=2*DRR_mu/mu_turb**2

solovevgeo=0
if solovevgeo==1:
    
    R0=6.2
    Z0=0

    R=np.linspace(4,8.5,120)
    Z=np.linspace(-4,4,200)

    r, z = np.meshgrid(R, Z)
    
    r=r/R0
    z=z/R0
    
    r2=r*r
    r3=r2*r
    r4=r3*r
    r5=r4*r
    r6=r5*r
    z2=z*z
    z3=z2*z
    z4=z3*z
    z5=z4*z
    z6=z5*z
    logr=np.log(r)
    
    B0=5.3
    psifac=200

    C = np.array([
             2.21808016e-02,  -1.28841781e-01,  -4.17718173e-02,
            -6.22680280e-02,   6.20083978e-03,  -1.20524711e-03,
            -3.70147050e-05,   0.00000000e+00,   0.00000000e+00,
             0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
            -0.155])
    
    psi = psifac * (
             (1-C[12]) * (r4/8)
              + C[12] * (r2*logr/2)
              + C[0]  * (1)
              + C[1]  * (r2)
              + C[2]  * (r2*logr - z2)
              + C[3]  * (r4 - 4*r2*z2)
              + C[4]  * (3*r4*logr - 9*r2*z2 - 12*r2*logr*z2 + 2*z4)
              + C[5]  * (r6 - 12*r4*z2 + 8*r2*z4)
              + C[6]  * (8*z6 - 140*r2*z4 - 120*r2*logr*z4 + 180*r4*logr*z2
                         + 75*r4*z2 - 15*r6*logr)
              + C[7]  * (z)
              + C[8]  * (z*r2)
              + C[9]  * (z3 - 3*z*r2*logr)
              + C[10] * (3*z*r4 - 4*z3*r2)
              + C[11] * (8*z5 - 45*z*r4 - 80*z3*r2*logr + 60*z*r4*logr) )
    
    BR =  -psifac/(r*R0*R0) * (
        C[2]  * (-2*z)
        + C[3]  * (-8*r2*z)
        + C[4]  * (-18*r2*z - 24*r2*logr*z + 8*z3)
        + C[5]  * (-24*r4*z + 32*r2*z3)
        + C[6]  * (48*z5 - 560*r2*z3 - 480*r2*logr*z3 +360*r4*logr*z
                   + 150*r4*z)
        + C[7]  * (1)
        + C[8]  * (r2)
        + C[9]  * (3*z2 - 3*r2*logr)
        + C[10] * (3*r4 - 12*z2*r2)
        + C[11] * (40*z4 - 45*r4 - 240*z2*r2*logr + 60*r4*logr) )
    
    BPHI=B0/r
    
    BZ =  psifac/(r*R0*R0) * (
        (1-C[12]) * (r3/2)
        + C[12] * (r/2 + r*logr)
        + C[1]  * (2*r)
        + C[2]  * (2*r*logr + r)
        + C[3]  * (4*r3 - 8*r*z2)
        + C[4]  * (12*r3*logr + 3*r3 - 30*r*z2 - 24*r*logr*z2)
        + C[5]  * (6*r5 - 48*r3*z2 + 16*r*z4)
        + C[6]  * (-400*r*z4 -240*r*logr*z4 + 720*r3*logr*z2 + 480*r3*z2
                   -90*r5*logr - 15*r5)
        + C[8]  * (2*z*r)
        + C[9]  * (-6*z*r*logr - 3*z*r)
        + C[10] * (12*z*r3 - 8*z3*r)
        + C[11] * (-120*z*r3-160*z3*r*logr-80*z3*r+240*z*r3*logr))

psi_lim=np.min(psi[:,0])
psi_min=np.min(psi[:])

new_paramaterization=1
if new_paramaterization==1:
    # A, mu[m], sigma[m]  —  model: sum A*exp(-0.5*((r-mu)/sigma)^2)
    gaussians = [
        (1.930982e-10, 0.06591262, 0.01758854),
        (8.701893e-09, 0.09884386, 0.01758854),
        (7.444911e-08, 0.13175743, 0.01758854),
        (4.521507e-08, 0.16459691, 0.01758854),
        (1.307995e-07, 0.22977723, 0.01758854),
        (5.733608e-07, 0.26193569, 0.01758854),
        (2.813727e-07, 0.29365476, 0.01758854),
        (2.960201e-07, 0.32480509, 0.01758854),
        (1.686966e-08, 0.35525330, 0.01758854),
        (8.048300e-07, 0.38487346, 0.01758854),
        (1.329798e-08, 0.44098735, 0.01758854),
        (5.040872e-07, 0.49172032, 0.01758854),
        (6.442105e-08, 0.53494040, 0.01758854),
    ]
    
    r     = np.linspace(0.05, 0.55, 500)
    total = sum(A*np.exp(-0.5*((r-mu)/s)**2) for A, mu, s in gaussians)
    
    plt.figure(figsize=(8,5))
    for A, mu, s in gaussians:
        plt.plot(r, A*np.exp(-0.5*((r-mu)/s)**2), '--', lw=0.8)  # individual components
    plt.plot(r, total, 'k-', lw=2, label='sum')
    plt.xlabel('minor radius r [m]')
    plt.ylabel(r'$\delta B_{norm}^2$')
    plt.legend()
    plt.tight_layout()
    plt.show()

#%% plotting

plt.rcParams['text.usetex'] = True
plt.rcParams['contour.negative_linestyle'] = 'solid'

SMALL_SIZE = 10
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=SMALL_SIZE)
plt.rc('figure', titlesize=SMALL_SIZE)

plot_psi=0
plot_b0=0
plot_rvskap=0
plot_br=0
plot_tmp=1

if plot_tmp==1:
    fig,ax=plt.subplots()
    
    #tmpfld=Lc.copy()
    #tmpfld=dBr_norm_squared.copy()
    #tmpfld=DRR.copy()
    #tmpfld=mu_drift.copy()
    #tmpfld=sigma.copy()
    #tmpfld=nu_drift.copy()
    tmpfld=nu_RR.copy()
    
    
    tmpfld[np.isnan(tmpfld)]=0

    mintmp=np.min(tmpfld)
    maxtmp=np.max(tmpfld)
    
    plot_bwr=0
    if plot_bwr==1:
        if maxtmp>np.abs(mintmp):
            mintmp=-maxtmp
        else:
            maxtmp=np.abs(mintmp)
        
    nlevs=50
    levs=np.linspace(mintmp,maxtmp,nlevs)
    ticklabels=np.linspace(mintmp,maxtmp,7)
    
    if plot_bwr==1:
        ct=ax.contourf(R,Z,tmpfld,levs,cmap='bwr')
    else:
        ct=ax.contourf(R,Z,tmpfld,levs)
        
    plt.plot(limR,limZ,color='k',linewidth=3)

    cbar=plt.colorbar(ct,ticks=ticklabels,format="%4.1e")
    #cbar.set_label('$B_r$', fontsize=12)

    plt.gca().set_aspect('equal')
    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    

    plt.savefig("tmp_plot.png", format="png", bbox_inches="tight")
    plt.show() 

if plot_rvskap==1:
    fig,ax=plt.subplots()
    
    plt.plot(radii,kappa_r)
    plt.plot(radii,popt[0]*(1+(radii/popt[1])**popt[2]))

if plot_b0==1:
    fig,ax=plt.subplots(1,3,figsize=(12,4))
    
    tmpfld=BR.copy()
    
    tmpfld[np.isnan(tmpfld)]=0

    mintmp=np.min(tmpfld)
    maxtmp=np.max(tmpfld)
    
    plot_bwr=1
    if plot_bwr==1:
        if maxtmp>np.abs(mintmp):
            mintmp=-maxtmp
        else:
            maxtmp=np.abs(mintmp)
        
    nlevs=50
    levs=np.linspace(mintmp,maxtmp,nlevs)
    ticklabels=np.linspace(mintmp,maxtmp,7)
    
    if plot_bwr==1:
        ct=ax[0].contour(R,Z,tmpfld,levs,cmap='bwr')
    else:
        ct=ax[0].contour(R,Z,tmpfld,levs)
        
    #plt.plot(limR,limZ,color='k',linewidth=3)
    ax[0].contour(R,Z,psi,levels=[psi_lim],colors='k',lw=2)

    cbar=plt.colorbar(ct,ax=ax[0],ticks=ticklabels,format="%4.1f")
    cbar.set_label('$B_R$', fontsize=12)

    ax[0].set_aspect('equal')
    ax[0].set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax[0].grid()

    tmpfld=BPHI.copy()
    
    tmpfld[np.isnan(tmpfld)]=0

    mintmp=np.min(tmpfld)
    maxtmp=np.max(tmpfld)
    
    plot_bwr=0
    if plot_bwr==1:
        if maxtmp>np.abs(mintmp):
            mintmp=-maxtmp
        else:
            maxtmp=np.abs(mintmp)
        
    nlevs=50
    levs=np.linspace(mintmp,maxtmp,nlevs)
    ticklabels=np.linspace(mintmp,maxtmp,7)
    
    if plot_bwr==1:
        ct=ax[1].contour(R,Z,tmpfld,levs,cmap='bwr')
    else:
        ct=ax[1].contour(R,Z,tmpfld,levs)
     
    ax[1].contour(R,Z,psi,levels=[psi_lim],colors='k',lw=2)
    #plt.plot(limR,limZ,color='k',linewidth=3)

    cbar=plt.colorbar(ct,ax=ax[1],ticks=ticklabels,format="%4.1f")
    cbar.set_label('$B_\\phi$', fontsize=12)

    ax[1].set_aspect('equal')
    ax[1].set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax[1].grid()
    
    tmpfld=BZ.copy()
    
    tmpfld[np.isnan(tmpfld)]=0

    mintmp=np.min(tmpfld)
    maxtmp=np.max(tmpfld)
    
    plot_bwr=1
    if plot_bwr==1:
        if maxtmp>np.abs(mintmp):
            mintmp=-maxtmp
        else:
            maxtmp=np.abs(mintmp)
        
    nlevs=50
    levs=np.linspace(mintmp,maxtmp,nlevs)
    ticklabels=np.linspace(mintmp,maxtmp,7)
    
    if plot_bwr==1:
        ct=ax[2].contour(R,Z,tmpfld,levs,cmap='bwr')
    else:
        ct=ax[2].contour(R,Z,tmpfld,levs)
        
    ax[2].contour(R,Z,psi,levels=[psi_lim],colors='k',lw=2)
        
    #plt.plot(limR,limZ,color='k',linewidth=3)

    cbar=plt.colorbar(ct,ax=ax[2],ticks=ticklabels,format="%4.1f")
    cbar.set_label('$B_Z$', fontsize=12)

    ax[2].set_aspect('equal')
    ax[2].set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax[2].grid()

    plt.savefig("B0_plot.pdf", format="pdf", bbox_inches="tight")
    plt.show() 

if plot_psi==1:
    fig,ax=plt.subplots()
    
    tmpfld=psi.copy()
    
    tmpfld[np.isnan(tmpfld)]=0

    mintmp=np.min(tmpfld)
    maxtmp=np.max(tmpfld)
    
    plot_bwr=0
    if plot_bwr==1:
        if maxtmp>np.abs(mintmp):
            mintmp=-maxtmp
        else:
            maxtmp=np.abs(mintmp)
        
    nlevs=50
    levs=np.linspace(mintmp,maxtmp,nlevs)
    ticklabels=np.linspace(mintmp,maxtmp,7)
    
    if plot_bwr==1:
        ct=ax.contour(R,Z,tmpfld,levs,cmap='bwr')
    else:
        ct=ax.contour(R,Z,tmpfld,levs)
        
    ax.contour(R,Z,psi,levels=[psi_lim],colors='k',lw=2)
        
    #plt.plot(limR,limZ,color='k',linewidth=3)

    cbar=plt.colorbar(ct,ticks=ticklabels,format="%4.1f")
    cbar.set_label('$\psi$', fontsize=12)

    plt.gca().set_aspect('equal')
    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    

    plt.savefig("psi_plot.pdf", format="pdf", bbox_inches="tight")
    plt.show() 
    
if plot_br==1:
    fig,ax=plt.subplots()
    
    tmpfld=Br.copy()**2
    
    tmpfld[np.isnan(tmpfld)]=0

    mintmp=np.min(tmpfld)
    maxtmp=np.max(tmpfld)
    
    plot_bwr=0
    if plot_bwr==1:
        if maxtmp>np.abs(mintmp):
            mintmp=-maxtmp
        else:
            maxtmp=np.abs(mintmp)
        
    nlevs=50
    levs=np.linspace(mintmp,maxtmp,nlevs)
    ticklabels=np.linspace(mintmp,maxtmp,7)
    
    if plot_bwr==1:
        ct=ax.contourf(R,Z,tmpfld,levs,cmap='bwr')
    else:
        ct=ax.contourf(R,Z,tmpfld,levs)
        
    plt.plot(limR,limZ,color='k',linewidth=3)

    cbar=plt.colorbar(ct,ticks=ticklabels,format="%4.1e")
    cbar.set_label('$B_r$', fontsize=12)

    plt.gca().set_aspect('equal')
    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    

    plt.savefig("br_plot.png", format="png", bbox_inches="tight")
    plt.show() 