import os
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import scipy as sci
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator

mu0=4*np.pi*10**(-7)
c = 2.99792458E8 #% Speed of light (m/s)
qe = 1.60217662E-19 #% Electron charge (C)
me = 9.10938356E-31 #% Electron mass (kg)
ep0 = 8.854E-12 #% Electric permitivity [C**2/(N*m**2)]
re = qe**2/(4*np.pi*ep0*me*c**2) #% classical electron radius (m)

#%% open Eric's matlab file for DIII-D magnetic field

on_nersc=1
if on_nersc==0:
    filename='/home/21b/ExperimentalInputs/DIIID/200236/Matts_Bfield_data.mat'

    mat=sci.io.loadmat(filename)

    Rlim=mat['RLim']
    Zlim=mat['ZLim']

    Rlim[3]=Rlim[2]
    Rlim[4]=Rlim[2]
    Rlim[5]=Rlim[2]
    Rlim[6]=Rlim[2]

#%%loading data from KORC output file

dir_num=1

#run_directory=['D3D_200236_MARS_CaseA_TEST18d1','D3D_200236_MARS_CaseA_TEST18d','D3D_200236_MARS_CaseA_TEST18d3','D3D_200236_MARS_CaseA_TEST18d4']
#run_directory=['../LOCAL/TEST8/OUT']
#run_directory=['GCeqn_GPU_TEST20b3']
#run_directory=['/home/21b/KORC_RUNS/FROM_PERLMUTTER/TEST20b3']
#run_directory=['../test/fio_m3dc1/tmp']
#run_directory=['/home/21b/KORC_RUNS/FROM_PERLMUTTER/TEST21']
#run_directory=['/home/21b/KORC_RUNS/FROM_PERLMUTTER/TEST20b_rr']
run_directory=['/pscratch/sd/m/mbeidler/KORC_GPU_RUNS/DIIID_177031_GPU_TEST21a']

for kk in range(0,dir_num):

    #filename=r"../"+run_directory[kk]+"/simulation_parameters.h5"
    filename=run_directory[kk]+"/simulation_parameters.h5"
    
    with h5py.File(filename,'r') as f:
        num_snapshots=f['simulation']['num_snapshots'][0]+1  
        nmpi=f['simulation']['nmpi'][0]
        output_cadence=f['simulation']['output_cadence'][0]
        t_steps_tot=f['simulation']['t_steps'][0]
        ppp=f['species']['ppp'][0]
        orbit_model=[x.decode() for x in f['simulation']['orbit_model']][0]
        field_eval=[x.decode() for x in f['simulation']['field_eval']][0]
        field_model=[x.decode() for x in f['simulation']['field_model']][0]
        outputs_list=[x.decode() for x in f['simulation']['outputs_list']][:]
        Bnorm=f['scales']['B'][0]
        tnorm=f['scales']['t'][0]
        Lnorm=f['scales']['l'][0]
        if field_model != 'M3D_C1':
            B0=f['fields']['Bo'][0]
            E0=f['fields']['Eo'][0]
            R0=f['fields']['Ro'][0]
            Z0=f['fields']['Zo'][0]
            if field_eval == 'eqn':
                ar=f['species']['ro'][0]
                a=f['fields']['a'][0]
                lam=f['fields']['lambda'][0]
                kappa=f['fields']['kappa'][0]
                q0=f['fields']['qo'][0]
            if field_eval != 'eqn':
                Rm=f['fields']['R'][:]
                Zm=f['fields']['Z'][:]
                #PSIPm=f['fields']['psi_p'][:]
                #BPHI1_Rem=f['fields']['BPHI1_Re'][:]
                
    if field_model == 'M3D_C1':
        filename=r"../"+run_directory[kk]+"/../m3dc1_outputs/C1.h5"
        with h5py.File(filename,'r') as f:
            Z0=f['scalars']['zmag'][0]
            R0=f['scalars']['xmag'][0]

        
    if num_snapshots<1:
        num_snapshots=1
    
    #filename=r"../"+run_directory[kk]+"/file_0.h5"
    filename=run_directory[kk]+"/file_0.h5"
    
    t_steps=np.zeros(num_snapshots,dtype=np.uint32) 
    timetmp=np.zeros((num_snapshots))

    with h5py.File(filename,'r') as f:
        for ii in range(0,num_snapshots):
            
            if kk==0:
                if ii==0:
                    t_steps[ii]=output_cadence*ii
                else:
                    if t_steps_tot<output_cadence:
                        t_steps[ii]=t_steps_tot
                    else:
                        t_steps[ii]=output_cadence*ii
                            
            else:
                t_steps[ii]=output_cadence*(ii+1)
                
            try:
                timetmp[ii]=f[str(t_steps[ii])]['time'][0]
            except:
                ii=ii-1
                break
            
    num_snapshots=ii+1
    if num_snapshots<1:
        num_snapshots=1
    
    nRE0=ppp*nmpi
    
    t_steps=np.zeros(num_snapshots,dtype=np.uint32)
    timetmp=np.zeros((num_snapshots))
    xxtmp=np.zeros((num_snapshots,nRE0))
    yytmp=np.zeros((num_snapshots,nRE0))
    Rtmp=np.zeros((num_snapshots,nRE0))
    PHItmp=np.zeros((num_snapshots,nRE0))
    zztmp=np.zeros((num_snapshots,nRE0))
    R0tmp=np.zeros((num_snapshots,nRE0))
    PHI0tmp=np.zeros((num_snapshots,nRE0))
    zz0tmp=np.zeros((num_snapshots,nRE0))
    vxtmp=np.zeros((num_snapshots,nRE0))
    vytmp=np.zeros((num_snapshots,nRE0))
    vztmp=np.zeros((num_snapshots,nRE0))
    pplltmp=np.zeros((num_snapshots,nRE0))
    mutmp=np.zeros((num_snapshots,nRE0))
    bRtmp=np.zeros((num_snapshots,nRE0))
    bPHItmp=np.zeros((num_snapshots,nRE0))
    bXtmp=np.zeros((num_snapshots,nRE0))
    bYtmp=np.zeros((num_snapshots,nRE0))
    bZtmp=np.zeros((num_snapshots,nRE0))
    curlbRtmp=np.zeros((num_snapshots,nRE0))
    curlbPHItmp=np.zeros((num_snapshots,nRE0))
    curlbZtmp=np.zeros((num_snapshots,nRE0))
    psiPtmp=np.zeros((num_snapshots,nRE0))
    eRtmp=np.zeros((num_snapshots,nRE0))
    ePHItmp=np.zeros((num_snapshots,nRE0))
    eZtmp=np.zeros((num_snapshots,nRE0))
    gtmp=np.zeros((num_snapshots,nRE0))
    etatmp=np.zeros((num_snapshots,nRE0))
    flagContmp=np.zeros((num_snapshots,nRE0))
    flagColtmp=np.zeros((num_snapshots,nRE0))
    flagREtmp=np.zeros((num_snapshots,nRE0))
    
    for jj in range(0,nmpi):
        
        #filename=r"../"+run_directory[kk]+"/file_"+str(jj)+".h5"
        filename=run_directory[kk]+"/file_"+str(jj)+".h5"
    
        with h5py.File(filename,'r') as f:
            for ii in range(0,num_snapshots):
                if jj==0:
                        
                    if kk==0:
                        if ii==0:
                            t_steps[ii]=output_cadence*ii
                        else:
                            if t_steps_tot<output_cadence:
                                t_steps[ii]=t_steps_tot
                            else:
                                t_steps[ii]=output_cadence*ii
                                    
                    else:
                        t_steps[ii]=output_cadence*(ii+1)
                    
                    timetmp[ii]=f[str(t_steps[ii])]['time'][0]
                    
                xxtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['X'][0][:]
                yytmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['X'][1][:]
                Rtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['Y'][0][:]
                PHItmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['Y'][1][:]
                if 'Y0' in outputs_list:
                    R0tmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['Y0'][0][:]
                    PHI0tmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['Y0'][1][:]
                bZtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['B'][2][:]
                if 'PSIp' in outputs_list:
                    psiPtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['PSIp'][:]
                if 'E' in outputs_list:
                    eRtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['E'][0][:]
                    ePHItmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['E'][1][:]
                    eZtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['E'][2][:]
                gtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['g'][:]
                etatmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['eta'][:]
                flagContmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['flagCon'][:]
                flagColtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['flagCol'][:]
                if 'curlb' in outputs_list:
                    curlbRtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['curlb'][0][:]
                    curlbPHItmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['curlb'][1][:]
                    curlbZtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['curlb'][2][:]
                if 'flagRE' in outputs_list:
                    flagREtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['flagRE'][:]
        
                if orbit_model=='FO':
                    zztmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['X'][2][:]
                    bXtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['B'][0][:]
                    bYtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['B'][1][:]
                    vxtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['V'][0][:]
                    vytmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['V'][1][:]
                    vztmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['V'][2][:]
                    
                elif orbit_model[0:2]=='GC':
                    if 'Y0' in outputs_list:
                        zz0tmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['Y0'][2][:]
                    zztmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['Y'][2][:]
                    bRtmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['B'][0][:]
                    bPHItmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['B'][1][:]
                    pplltmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['V'][0][:]
                    mutmp[ii,jj*ppp:(jj+1)*ppp]=f[str(t_steps[ii])]['spp_1']['V'][1][:]
                    
    if kk==0:
        time=timetmp
        xx=xxtmp
        yy=yytmp
        if 'Y0' in outputs_list:
            R00=R0tmp
            PHI00=PHI0tmp
            Z00=zz0tmp
        R=Rtmp
        PHI=PHItmp
        zz=zztmp
        bZ=bZtmp
        if 'PSIp' in outputs_list: 
            psiP=psiPtmp
        eR=eRtmp
        ePHI=ePHItmp
        eZ=eZtmp
        g=gtmp
        eta=etatmp
        flagCon=flagContmp
        flagCol=flagColtmp
        if 'curlb' in outputs_list:
            curlbR=curlbRtmp
            curlbPHI=curlbPHItmp
            curlbZ=curlbZtmp
        if 'flagRE' in outputs_list:
            flagRE=flagREtmp
        if orbit_model=='FO':
            bX=bXtmp
            bY=bYtmp
            vx=vxtmp
            vy=vytmp
            vz=vztmp
            
        elif orbit_model[0:2]=='GC':
            bR=bRtmp
            bPHI=bPHItmp
            ppll=pplltmp
            mu=mutmp
    else:
        
        add_time=1
        add_part=0
        
        if add_time==1:

            time=np.concatenate((time,timetmp-time[-1]))
            xx=np.concatenate((xx,xxtmp),axis=0)
            yy=np.concatenate((yy,yytmp),axis=0)
            if 'Y0' in outputs_list:
                R00=np.concatenate((R00,R0tmp),axis=0)
                PHI0=np.concatenate((PHI00,PHI0tmp),axis=0)
                Z00=np.concatenate((Z00,zz0tmp),axis=0)
            R=np.concatenate((R,Rtmp),axis=0)
            PHI=np.concatenate((PHI,PHItmp),axis=0)
            zz=np.concatenate((zz,zztmp),axis=0)
            bZ=np.concatenate((bZ,bZtmp),axis=0)
            if 'PSIp' in outputs_list:
                psiP=np.concatenate((psiP,psiPtmp),axis=0)
            eR=np.concatenate((eR,eRtmp),axis=0)
            ePHI=np.concatenate((ePHI,ePHItmp),axis=0)
            eZ=np.concatenate((eZ,eZtmp),axis=0)
            g=np.concatenate((g,gtmp),axis=0)
            eta=np.concatenate((eta,etatmp),axis=0)
            flagCon=np.concatenate((flagCon,flagContmp),axis=0)
            flagCol=np.concatenate((flagCol,flagColtmp),axis=0)
            if 'curlb' in outputs_list:
                curlbR=np.concatenate((curlbR,curlbRtmp),axis=0)
                curlbPHI=np.concatenate((curlbPHI,curlbPHItmp),axis=0)
                curlbZ=np.concatenate((curlbZ,curlbZtmp),axis=0)
            if 'flagRE' in outputs_list:
                flagRE=np.concatenate((flagRE,flagREtmp),axis=0)
                
            if orbit_model=='FO':
                bX=np.concatenate((bX,bXtmp),axis=0)
                bY=np.concatenate((bY,bYtmp),axis=0)
                vx=np.concatenate((vx,vxtmp),axis=0)
                vy=np.concatenate((vy,vytmp),axis=0)
                vz=np.concatenate((vz,vztmp),axis=0)
                
            elif orbit_model[0:2]=='GC':
                bR=np.concatenate((bR,bRtmp),axis=0)
                bPHI=np.concatenate((bPHI,bPHItmp),axis=0)
                ppll=np.concatenate((ppll,pplltmp),axis=0)
                mu=np.concatenate((mu,mutmp),axis=0)
        
        if add_part==1:
            min_steps = min(time.shape, timetmp.shape)[0]
            
            time=time[:min_steps]
            xx=np.concatenate((xx[:min_steps, :],xxtmp[:min_steps, :]),axis=1)
            yy=np.concatenate((yy[:min_steps, :],yytmp[:min_steps, :]),axis=1)
            R=np.concatenate((R[:min_steps, :],Rtmp[:min_steps, :]),axis=1)
            PHI=np.concatenate((PHI[:min_steps, :],PHItmp[:min_steps, :]),axis=1)
            zz=np.concatenate((zz[:min_steps, :],zztmp[:min_steps, :]),axis=1)
            bZ=np.concatenate((bZ[:min_steps, :],bZtmp[:min_steps, :]),axis=1)
            if 'PSIp' in outputs_list:
                psiP=np.concatenate((psiP[:min_steps, :],psiPtmp[:min_steps, :]),axis=1)
            eR=np.concatenate((eR[:min_steps, :],eRtmp[:min_steps, :]),axis=1)
            ePHI=np.concatenate((ePHI[:min_steps, :],ePHItmp[:min_steps, :]),axis=1)
            eZ=np.concatenate((eZ[:min_steps, :],eZtmp[:min_steps, :]),axis=1)
            g=np.concatenate((g[:min_steps, :],gtmp[:min_steps, :]),axis=1)
            eta=np.concatenate((eta[:min_steps, :],etatmp[:min_steps, :]),axis=1)
            flagCon=np.concatenate((flagCon[:min_steps, :],flagContmp[:min_steps, :]),axis=1)
            flagCol=np.concatenate((flagCol[:min_steps, :],flagColtmp[:min_steps, :]),axis=1)
            
            if orbit_model=='FO':
                bX=np.concatenate((bX[:min_steps, :],bXtmp[:min_steps, :]),axis=1)
                bY=np.concatenate((bY[:min_steps, :],bYtmp[:min_steps, :]),axis=1)
                vx=np.concatenate((vx[:min_steps, :],vxtmp[:min_steps, :]),axis=1)
                vy=np.concatenate((vy[:min_steps, :],vytmp[:min_steps, :]),axis=1)
                vz=np.concatenate((vz[:min_steps, :],vztmp[:min_steps, :]),axis=1)
                
            elif orbit_model[0:2]=='GC':
                bR=np.concatenate((bR[:min_steps, :],bRtmp[:min_steps, :]),axis=1)
                bPHI=np.concatenate((bPHI[:min_steps, :],bPHItmp[:min_steps, :]),axis=1)
                ppll=np.concatenate((ppll[:min_steps, :],pplltmp[:min_steps, :]),axis=1)
                mu=np.concatenate((mu[:min_steps, :],mutmp[:min_steps, :]),axis=1)

#%% open KORC field file

useextfield=0

if useextfield==1:

    #filename=r"/home/21b/KORC_RUNS/LOCAL/TEST2/AORSA_D3D_171089_200MHz_EFITgrid.h5"
    
    if on_nersc==0:
        filename=r"../"+run_directory[0]+"/DIII_200236_MARS_CaseA.h5"
    else:
        filename=r"../../"+run_directory[0]+"/DIII_200236_MARS_CaseA.h5"
    
    with h5py.File(filename,'r') as f:
        NR=int(f['NR'][0])           
        NZ=int(f['NZ'][0])
        
        R0=f['Ro'][0]     
        Z0=f['Zo'][0]      
        B0=f['Bo'][0]
        
        Rg=f['R'][:]
        Zg=f['Z'][:]
        
        Rg=np.array(Rg)
        Zg=np.array(Zg)
        
        PSIp=f['PSIp'][:]
        FLAG=f['FLAG'][:]
        
        PSIp=np.array(PSIp)
        FLAG=np.array(FLAG)
    
    
    ZZg,RRg=np.meshgrid(Zg,Rg,indexing='ij')
    
    timeind=39
    
    interp_psip = RegularGridInterpolator((Rg,Zg), np.transpose(PSIp),method='cubic')
    
    dpsipdz,dpsipdr=np.gradient(PSIp,Zg,Rg)
    
    BR_g=-dpsipdz*np.divide(1,RRg)
    BZ_g=dpsipdr*np.divide(1,RRg)
    BPHI_g=-B0*R0/RRg#np.divide(1,RRg)
    
    interp_br = RegularGridInterpolator((Rg,Zg), np.transpose(BR_g[:,:]),method='cubic')
    interp_bz = RegularGridInterpolator((Rg,Zg), np.transpose(BZ_g[:,:]),method='cubic')
    
    Rq=1.7
    Zq=0.
    
    psip_q=interp_psip([Rq,Zq])
    br_q=interp_br([Rq,Zq])
    bz_q=interp_bz([Rq,Zq])

#%% open Eric's matlab file for RE kinetic distributions

usekindist=0

if usekindist==1:
    
    if on_nersc==0:
        filename=r"../"+run_directory[0]+"/Hollmann_PDF_DIII_200236.h5"
    else:
        filename=r"../../"+run_directory[0]+"/Hollmann_PDF_DIII_200236.h5"
    
    with h5py.File(filename,'r') as f:
        EfitVc=f['E'][:]/10**6
        rhoVRE=f['rho'][:]
        fM=np.transpose(f['fRE_E'][:])
        #costhA=mat['saveA']['fRE_pitch'][0]
        #densREV=10**19*mat['saveA']['densREV'][0][0][0]
    
    #filename='/home/21b/ExperimentalInputs/DIIID/200236/Matts_fM_data.mat'
    
    #mat=sci.io.loadmat(filename)
    
    #EfitVc=mat['saveA']['EfitVc'][0][0][0]
    #rhoVRE=mat['saveA']['rhoVRE'][0][0][0]
    #rhoV=mat['saveA']['rhoV'][0][0][0]
    #fM=mat['saveA']['fM'][0][0]
    #costhA=mat['saveA']['costhA'][0][0]
    #densREV=10**19*mat['saveA']['densREV'][0][0][0]
    
    #NE=np.shape(EfitVc)[0]
    #Nrho=np.shape(rhoVRE)[0]
    
    fMnorm=np.trapz(fM,EfitVc*10**6)

#%% Analysis

if orbit_model=='FO':
    R=np.sqrt(xx**2+yy**2)
    PHI=np.arctan2(yy,xx)
    
    bR=bX*np.cos(PHI)+bY*np.sin(PHI)
    bPHI=bY*np.cos(PHI)-bX*np.sin(PHI)
    
    vmag=np.sqrt(vx**2+vy**2+vz**2)
    bmag=np.sqrt(bX**2+bY**2+bZ**2)
elif orbit_model[0:2]=='GC':
    
    bX=bR*np.cos(PHI)-bPHI*np.sin(PHI)
    bY=bR*np.sin(PHI)+bPHI*np.cos(PHI)
    
    bmag=np.sqrt(bR**2+bPHI**2+bZ**2)
    bpol=np.sqrt(bR**2+bZ**2)
    vmag=c*np.sqrt(1-(1/g)**2)
    pmag=me*c*np.sqrt(g**2-1)
    vperp=vmag*np.sin(np.radians(eta))
    vpll=vmag*np.cos(np.radians(eta))

rm=np.sqrt((R-R0)**2+(zz-Z0)**2)

flagDecon=np.zeros(np.shape(flagCon))
flagDecon[(flagCon<1) & (flagRE>0)]=1
confined=np.sum(flagCon,axis=1)
deconfined=np.sum(flagDecon,axis=1)

flagTherm=np.zeros(np.shape(flagCol))
flagTherm[flagCol<1 & (flagRE>0)]=1
Thermal=np.sum(flagTherm,axis=1)
Energetic=np.sum(flagCol,axis=1)

if 'flagRE' in outputs_list:
  Total=np.sum(flagRE,axis=1)
    
  flagActive=flagRE*flagCon*flagCol
  flagActive[np.where(np.isnan(R))]=0
  Active=np.sum(flagActive,axis=1)

  flagPrimary=np.zeros(np.shape(flagRE))
  flagPrimary[:,flagRE[0,:]>0]=1

  flagSecondary=np.zeros(np.shape(flagRE))
  flagSecondary[:,flagRE[0,:]<1]=1
  flagSecondary[:,flagRE[-1,:]<1]=0

Ipart=qe*vpll*bPHI/(2*np.pi*R*bmag)
Ipart[flagActive==0]=0
Itot=np.sum(Ipart,axis=1)
Ipri=np.sum(Ipart*flagPrimary,axis=1)
Isec=np.sum(Ipart*flagSecondary,axis=1)

IPOLpart=qe*vpll*bpol/(2*np.pi*R*bmag)
IPOLpart[flagActive==0]=0
IPOL=np.sum(IPOLpart,axis=1)


#KE=me*c**2/qe*np.sum((g-1)*flagActive)
#KEtot=me*c**2/qe*np.sum((1-1)*flagRE)

K=(g-1)*(me*c**2/qe)

GR=vperp/(qe*bmag/(me*g))

#%% analytic fields

if (field_model != 'M3D_C1') and (field_eval == 'eqn'):

    Rm=np.linspace(R0-a,R0+a,100)
    Zm=np.linspace(Z0-a,Z0+a,100)
    RRm, ZZm = np.meshgrid(Rm, Zm)
    
    qa=q0*(1+(a/lam)**2)
    rm=np.sqrt((RRm-R0)**2+(ZZm-Z0)**2)
    
    rm_elong=np.sqrt((RRm-R0)**2+(ZZm-Z0)**2/kappa**2)
    
    theta=np.atan2((ZZm-Z0),(RRm-R0))
    psi=lam**2*B0/(2*q0)*np.log(1+(rm_elong/lam)**2)
    
    limR=R0+a*np.cos(np.linspace(0,2*np.pi,100))/kappa
    limZ=Z0+a*np.sin(np.linspace(0,2*np.pi,100))

#%% Calculate and save facetted wall impacts

# 1. Base Configuration
nt = 42
Rc = 1.016
seg = 5000
        
nsam_chi=10
num_timesteps,num_particles = mu.shape
chi = np.zeros((num_timesteps,num_particles,nsam_chi))

# Regular polygon wall definition
philim = np.linspace(0, 2 * np.pi, seg)
rlim = Rc * np.cos(np.pi / nt) / (np.cos(np.mod(philim, 2 * np.pi / nt) - np.pi / nt))
xiwl = rlim * np.cos(philim)
yiwl = rlim * np.sin(philim)

plotinter = 0
if plotinter == 1:
    fig, ax = plt.subplots()
    ax.plot(xiwl, yiwl)
    ax.set_aspect('equal')

# 2. Initialize output arrays matching your exact matrix profiles
Rint = np.zeros_like(R)
PHIint = np.zeros_like(R)
Zint = np.zeros_like(R)
inc = np.zeros_like(R)
conlossind = np.zeros_like(flagDecon)
PFCinter = np.zeros_like(flagDecon)

# 1. Compute the time-step difference mask (yields 13 rows)
standard_loss = (flagDecon[1:, :] - flagDecon[:-1, :]) == 1
is_new_born = (flagRE[1:, :] - flagRE[:-1, :]) == 1
born_deconfined = is_new_born & (flagDecon[1:, :] == 1)
mask_13 = standard_loss | born_deconfined

# 2. Expand mask to 14 indices by prepending a row of False for the initial timestep
mask = np.zeros((num_snapshots, num_particles), dtype=bool)
mask[1:, :] = mask_13

# If any elements match the condition, compute the geometric vector slices instantly
if np.any(mask):
    
    conlossind[mask] = 1

    # Extract dynamic active element vectors
    r00_m   = R00[1:, :][mask[1:, :]]
    phi00_m = PHI00[1:, :][mask[1:, :]]
    r_m     = R[1:, :][mask[1:, :]]
    phi_m   = PHI[1:, :][mask[1:, :]]
    z00_m   = Z00[1:, :][mask[1:, :]]
    z_m     = zz[1:, :][mask[1:, :]]
    
    # Convert cylindrical paths to Cartesian coordinates
    x1 = r00_m * np.cos(phi00_m)
    y1 = r00_m * np.sin(phi00_m)
    x2 = r_m * np.cos(phi_m)
    y2 = r_m * np.sin(phi_m)

    # Wall segment angle tracking
    dphi = 2 * np.pi / nt
    mint = np.floor(phi00_m / dphi)
    minPHIt = mint * dphi
    maxPHIt = (mint + 1) * dphi

    # Standard normal vector vectors
    nhatx = np.cos((maxPHIt + minPHIt) / 2)
    nhaty = np.sin((maxPHIt + minPHIt) / 2)

    # Target polygon boundary vertices
    xt1 = Rc * np.cos(maxPHIt)
    yt1 = Rc * np.sin(maxPHIt)
    xt2 = Rc * np.cos(minPHIt)
    yt2 = Rc * np.sin(minPHIt)

    # Line intersection math determinants
    dx = x1 - x2
    dxt = xt1 - xt2
    dy = y1 - y2
    dyt = yt1 - yt2
    dr = np.sqrt(dx**2 + dy**2)

    D = x1 * y2 - y1 * x2
    Dt = xt1 * yt2 - yt1 * xt2
    Dint = dx * dyt - dy * dxt

    # Ray intersection coordinate maps
    xint = (D * dxt - Dt * dx) / Dint
    yint = (D * dyt - Dt * dy) / Dint

    # Calculate intersection bounds check
    is_outside = ((xint > x1) & (xint > x2)) | ((xint < x1) & (xint < x2))

    PFCinter[mask] = np.where(is_outside, 0, 1)

    # Filter out alternative structures (other wall losses)
    other_wall_loss = (PFCinter[mask] == 0)
    xint[other_wall_loss] = np.nan
    yint[other_wall_loss] = np.nan

    # Write metrics directly into main arrays (No manual column-slicing needed at the end!)
    Rint[mask]   = np.sqrt(xint**2 + yint**2)
    PHIint[mask] = np.arctan2(yint, xint)

    del_val = np.where(dx != 0, (xint - x1) / dx, 0)
    Zint[mask] = z00_m + del_val * (z00_m - z_m)
    Zint[mask][other_wall_loss] = np.nan  

    inc[mask] = (nhatx * dx / dr + nhaty * dy / dr)
    inc[mask][other_wall_loss] = np.nan

    # --- Debug checks corresponding to your loops ---
    # Obtain exact row indexes (j) where incidence is negative to replicate your tracking conditions
    rows, cols = np.where((inc < 0) & (conlossind == 1))
    for j, i in zip(rows, cols):
        REnum = np.mod(j, 35000)
        # Add your custom print or tracking outputs here if debugging...

    # --- Optional Diagnostic Plotting ---
    if plotinter == 1:
        for x_start, x_end, y_start, y_end, hit in zip(x1, x2, y1, y2, PFCinter[mask]):
            if hit ==1:
                ax.plot([x_start, x_end], [y_start, y_end], 'o-')
        plt.show()

# --- 4. DATA SELECTION & INCIDENT INDEX IDENTIFICATION ---
# Filter particles: Must have a loss event (conlossind == 1) AND an active wall intersection (inc > 0)
# FIX: axis=0 looks across all timesteps (rows) for each independent particle (columns)
valid_particles = np.any(conlossind == 1, axis=0) & np.any(inc > 0, axis=0)

# Isolate column (particle) indices passing this condition
# valid_particles is 1D with shape (num_particles,)
jj_idx = np.where(valid_particles)[0]

if len(jj_idx) > 0:
    # Find the FIRST row timestep index where the loss event occurs for each valid particle
    # FIX: axis=0 searches down the rows (timesteps) for the first True element
    indt_idx = np.argmax(conlossind[:, jj_idx] == 1, axis=0)
    N_incidents = len(jj_idx)
else:
    N_incidents = 0

# --- 5. VECTORIZED GYROPHASE (CHI) CALCULATION FOR INCIDENTS ---
# Initialize the final 3D chi framework with zeros
# CORRECTED CONFIGURATION SHAPE: (num_timesteps, num_particles, nsam_chi)
chi = np.zeros((num_timesteps, num_particles, nsam_chi))

if N_incidents > 0:
    # Extract only the physical coordinates at the exact timestamp of impact
    # FIX: Row index is indt_idx (timesteps), Column index is jj_idx (particles)
    inc_inc = inc[indt_idx, jj_idx]
    eta_inc = np.radians(eta[indt_idx, jj_idx])
    
        # Look up the flag values at the exact (timestep, particle) coordinate of impact
    # This yields a 1D array of shape (N_incidents,)
    is_primary_at_impact = flagPrimary[indt_idx, jj_idx]
    is_secondary_at_impact = flagSecondary[indt_idx, jj_idx]
    
    # Convert them to strict boolean masks
    # (Evaluates to True if the flag is 1 or True at the moment of impact)
    primary_mask = (is_primary_at_impact == 1)
    secondary_mask = (is_secondary_at_impact == 1)

    # Safeguard: Clip parameters to [-1, 1] bounds to prevent numerical arcsin NaN breaks
    inc_inc = np.clip(inc_inc, -1.0, 1.0)

    # 1D mathematical profiles
    term1 = -np.cos(np.arcsin(inc_inc)) * np.sin(eta_inc)
    term2 = inc_inc * np.abs(np.cos(eta_inc))

    # Replicate into a 2D grid layout: Shape (N_incidents, nsam_chi)
    term1_2d = np.repeat(term1[:, np.newaxis], nsam_chi, axis=1)
    term2_2d = np.repeat(term2[:, np.newaxis], nsam_chi, axis=1)

    # Initialize random sampling distribution
    chi_working = 2 * np.pi * np.random.rand(N_incidents, nsam_chi)
    incangle = term1_2d * np.sin(chi_working) + term2_2d
    invalid_mask = (incangle < 0) | (incangle > 1)

    # 2D Vectorized Rejection Sampling Loop
    while np.any(invalid_mask):
        num_needed = np.sum(invalid_mask)
        new_samples = 2 * np.pi * np.random.rand(num_needed)
        chi_working[invalid_mask] = new_samples
        
        # Recalculate only the failed validation checkpoints
        incangle[invalid_mask] = term1_2d[invalid_mask] * np.sin(new_samples) + term2_2d[invalid_mask]
        invalid_mask = (incangle < 0) | (incangle > 1)

    # Drop verified 2D samples back to the master 3D framework instantly
    # CORRECTED 3D ASSIGNMENT INDEXING: chi[row_indices, col_indices, :]
    chi[indt_idx, jj_idx, :] = chi_working
    
# --- 1. CONFIGURATION SETUP ---
collapse = 1  # Set to 1 to fold all particles onto a single wall segment

half_tile = np.pi / nt

# Setup bin parameters based on collapse toggle
if collapse == 0:
    nPHIbins = 200
    PHIbinedges = np.linspace(0, 2 * np.pi, nPHIbins + 1)
else:
    # Collapse folds everything onto a single tile domain: [0, 2*pi/nt]
    nPHIbins = 50
    #PHIbinedges = np.linspace(0, 2 * np.pi / nt, nPHIbins + 1)
    PHIbinedges = np.linspace(-half_tile, half_tile, nPHIbins + 1)

edge = np.linspace(0, 2 * np.pi, nt + 1)
dPHI = PHIbinedges[1] - PHIbinedges[0]

# Calculate bin centers directly without a loop
# (Equivalent to your MATLAB for i=1:size(PHIvals,2) loop)
PHIvals = (PHIbinedges[:-1] + PHIbinedges[1:]) / 2


# --- 2. VECTORIZED PHI BINDING FOR INCIDENT PARTICLES ---
# Using the indt_idx (timesteps) and jj_idx (particles) from your previous gyrophase block:
if N_incidents > 0:
    # Pull the calculated intersection angles for only the incident group
    # Shape: (N_incidents,)
    phi_incidents_raw = PHIint[indt_idx, jj_idx]

    half_tile = np.pi / nt
    dphi_face = 2 * np.pi / nt
    phi_wrapped = np.mod(phi_incidents_raw + half_tile, dphi_face) - half_tile

    Bx_raw = bX[indt_idx, jj_idx]
    By_raw = bY[indt_idx, jj_idx]
    Bz_raw = bZ[indt_idx, jj_idx]

    # Calculate native cylindrical magnetic components (Invariants across all tile faces)
    cos_raw_1d = np.cos(phi_incidents_raw)
    sin_raw_1d = np.sin(phi_incidents_raw)
    Br_field_local   =  Bx_raw * cos_raw_1d + By_raw * sin_raw_1d
    Bphi_field_local = -Bx_raw * sin_raw_1d + By_raw * cos_raw_1d

    # Project fields into the reference Cartesian frame of Tile 0
    cos_wrap_1d = np.cos(phi_wrapped)
    sin_wrap_1d = np.sin(phi_wrapped)
    Bx_collapsed = Br_field_local * cos_wrap_1d - Bphi_field_local * sin_wrap_1d
    By_collapsed = Br_field_local * sin_wrap_1d + Bphi_field_local * cos_wrap_1d


    # --- 3. REPLICATE AND FLATTEN ALL ARRAYS TO MATCH CHI (Shape: N_incidents * nsam_chi) ---
    # Extract raw parameters
    vmag_inc_1d = vmag[indt_idx, jj_idx]
    eta_inc_1d  = eta[indt_idx, jj_idx]
    chi_inc_flat = chi[indt_idx, jj_idx, :].ravel()  # Master dimension footprint

    # Repeat spatial angle profiles
    phi_raw_flat     = np.repeat(phi_incidents_raw[:, np.newaxis], nsam_chi, axis=1).ravel()
    phi_wrapped_flat = np.repeat(phi_wrapped[:, np.newaxis], nsam_chi, axis=1).ravel()

    # Repeat physical magnitudes and field arrays
    vmag_inc_flat = np.repeat(vmag_inc_1d[:, np.newaxis], nsam_chi, axis=1).ravel()
    eta_inc_flat  = np.repeat(eta_inc_1d[:, np.newaxis], nsam_chi, axis=1).ravel()
    
    Bx_raw_flat   = np.repeat(Bx_raw[:, np.newaxis], nsam_chi, axis=1).ravel()
    By_raw_flat   = np.repeat(By_raw[:, np.newaxis], nsam_chi, axis=1).ravel()
    Bz_raw_flat   = np.repeat(Bz_raw[:, np.newaxis], nsam_chi, axis=1).ravel()
    
    
def cylindrical_to_cartesian(r, phi, z):
    """
    Transforms 3D cylindrical coordinates (r, phi, z) to Cartesian (x, y, z).
    Accepts scalars, 1D arrays, 2D matrices, or 3D grids.
    """
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    # Z remains unchanged in cylindrical space
    return x, y, z

def velocity_magnetic_to_collapsed_cartesian(v, eta, chi, Bx, By, Bz, phi_raw, phi_wrap):
    """
    Transforms spherical velocities relative to B directly into a wall-collapsed
    Cartesian frame using local cylindrical transformations to isolate true transport vectors.
    """
    # Reconstruct total magnetic magnitude profile
    B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)
    B_mag_safe = np.where(B_mag == 0, 1e-10, B_mag)
    
    # Parallel tracking vector unit projections
    bx = Bx / B_mag_safe
    by = By / B_mag_safe
    bz = Bz / B_mag_safe
    
    # Build local radial unit vector at the particle's RAW impact position
    er_x = np.cos(phi_raw)
    er_y = np.sin(phi_raw)
    er_z = np.zeros_like(phi_raw)
    
    # Isolate tracking vector 1 perpendicular to the field line using the normal vector
    dot_b_er = bx * er_x + by * er_y + bz * er_z
    perp1_x = er_x - dot_b_er * bx
    perp1_y = er_y - dot_b_er * by
    perp1_z = er_z - dot_b_er * bz
    perp1_norm = np.sqrt(perp1_x**2 + perp1_y**2 + perp1_z**2)
    perp1_norm_safe = np.where(perp1_norm == 0, 1.0, perp1_norm)
    
    e2_x = perp1_x / perp1_norm_safe
    e2_y = perp1_y / perp1_norm_safe
    e2_z = perp1_z / perp1_norm_safe
    
    # Tracking vector 2 completes the orthogonal coordinate system (b x e2)
    e3_x = by * e2_z - bz * e2_y
    e3_y = bz * e2_x - bx * e2_z
    e3_z = bx * e2_y - by * e2_x
    
    # Project local spherical components into the customized coordinate system
    eta_rad = np.radians(eta)
    v_parallel = v * np.cos(eta_rad)
    v_perp1    = v * np.sin(eta_rad) * np.cos(chi)
    v_perp2    = v * np.sin(eta_rad) * np.sin(chi)
    
    # Reconstruct native un-collapsed Cartesian velocities
    vx_raw = v_parallel * bx + v_perp1 * e2_x + v_perp2 * e3_x
    vy_raw = v_parallel * by + v_perp1 * e2_y + v_perp2 * e3_y
    vz_out = v_parallel * bz + v_perp1 * e2_z + v_perp2 * e3_z
    
    # Convert raw Cartesian velocities to local Cylindrical coordinates (Invariants)
    cos_raw = np.cos(phi_raw)
    sin_raw = np.sin(phi_raw)
    vr_local   =  vx_raw * cos_raw + vy_raw * sin_raw
    vphi_local = -vx_raw * sin_raw + vy_raw * cos_raw
    
    # Project components into the Cartesian coordinates of your reference tile (phi_wrap)
    cos_wrap = np.cos(phi_wrap)
    sin_wrap = np.sin(phi_wrap)
    vx_collapsed = vr_local * cos_wrap - vphi_local * sin_wrap
    vy_collapsed = vr_local * sin_wrap + vphi_local * cos_wrap
    
    return vx_collapsed, vy_collapsed, vz_out

time_inc_1d = time[indt_idx]

inc_inc_1d  = inc[indt_idx, jj_idx]

Rint_inc_1d = Rint[indt_idx,jj_idx]
PHIint_inc_1d = phi_wrapped
Zint_inc_1d = Zint[indt_idx,jj_idx]

#bX_inc_1d  = bX[indt_idx, jj_idx]
#bY_inc_1d  = bY[indt_idx, jj_idx]

bX_inc_1d  = Bx_collapsed
bY_inc_1d  = By_collapsed

bR_inc_1d  = bR[indt_idx, jj_idx]
bPHI_inc_1d  = bPHI[indt_idx, jj_idx]
bZ_inc_1d  = bZ[indt_idx, jj_idx]

# Extract the 2D chi values (Shape: (N_incidents, nsam_chi))
chi_inc_2d  = chi[indt_idx, jj_idx, :].ravel()

# Replicate all 1D metrics into identical 2D grids (Shape: (N_incidents, nsam_chi))
time_inc_2d = np.repeat(time_inc_1d[:, np.newaxis], nsam_chi, axis=1).ravel()

vmag_inc_2d = np.repeat(vmag_inc_1d[:, np.newaxis], nsam_chi, axis=1).ravel()
eta_inc_2d  = np.repeat(eta_inc_1d[:, np.newaxis], nsam_chi, axis=1).ravel()
inc_inc_2d  = np.repeat(inc_inc_1d[:, np.newaxis], nsam_chi, axis=1).ravel()

Rint_inc_2d    = np.repeat(Rint_inc_1d[:, np.newaxis], nsam_chi, axis=1).ravel()
PHIint_inc_2d  = np.repeat(PHIint_inc_1d[:, np.newaxis], nsam_chi, axis=1).ravel()
Zint_inc_2d    = np.repeat(Zint_inc_1d[:, np.newaxis], nsam_chi, axis=1).ravel()

bX_inc_2d =  np.repeat(bX_inc_1d[:, np.newaxis], nsam_chi, axis=1).ravel()
bY_inc_2d  = np.repeat(bY_inc_1d[:, np.newaxis], nsam_chi, axis=1).ravel()
bR_inc_2d =  np.repeat(bR_inc_1d[:, np.newaxis], nsam_chi, axis=1).ravel()
bPHI_inc_2d  = np.repeat(bPHI_inc_1d[:, np.newaxis], nsam_chi, axis=1).ravel()
bZ_inc_2d  = np.repeat(bZ_inc_1d[:, np.newaxis], nsam_chi, axis=1).ravel()

xint, yint, zint = cylindrical_to_cartesian(Rint_inc_2d,PHIint_inc_2d,Zint_inc_2d)
xint_1d, yint_1d, zint_1d = cylindrical_to_cartesian(Rint_inc_1d,PHIint_inc_1d,Zint_inc_1d)

vx, vy, vz = velocity_magnetic_to_collapsed_cartesian(
        vmag_inc_flat, 
        eta_inc_flat, 
        chi_inc_flat+3*np.pi/2, 
        Bx_raw_flat, 
        By_raw_flat, 
        Bz_raw_flat, 
        phi_raw_flat, 
        phi_wrapped_flat
    )

sav=0
if sav==1:
       
    filename='IWL_impacts_TEST21.h5'

    with h5py.File(filename, "w") as f:
        dset=f.create_dataset('NRE',(1,),dtype='i')
        dset[0]=N_incidents*nsam_chi
        dset=f.create_dataset('Time',(N_incidents*nsam_chi,),dtype='f')
        dset[:]=time_inc_2d
        #dset=f.create_dataset('R',(N_incidents*nsam_chi,),dtype='f')
        #dset[:]=Rint_inc_2d
        #dset=f.create_dataset('PHI',(N_incidents*nsam_chi,),dtype='f')
        #dset[:]=PHIint_inc_2d
        dset=f.create_dataset('X',(N_incidents*nsam_chi,),dtype='f')
        dset[:]=xint
        dset=f.create_dataset('Y',(N_incidents*nsam_chi,),dtype='f')
        dset[:]=yint
        dset=f.create_dataset('Z',(N_incidents*nsam_chi,),dtype='f')
        dset[:]=zint
        dset=f.create_dataset('VX',(N_incidents*nsam_chi,),dtype='f')
        dset[:]=vx
        dset=f.create_dataset('VY',(N_incidents*nsam_chi,),dtype='f')
        dset[:]=vy
        dset=f.create_dataset('VZ',(N_incidents*nsam_chi),dtype='f')
        dset[:]=vz

#%% Save DiMES impacts

DiMESdepo=0
if DiMESdepo==1:

    theta_geo=np.arctan2(zz-Z0,R-R0)
    
    #DiMES impact calculation
    DiMESloc_cyl=[1.485,np.deg2rad(150),-1.245] #In (R,PHI,Z)
    
    theta_DiMES=np.arctan2(DiMESloc_cyl[2]-Z0,DiMESloc_cyl[0]-R0)
    
    #DiMESdims=[0.025,0.01] # (radius,height of dome) for semi-spheroid
    DiMESdims=[0.033,0.01] # (radius,height of dome) for section of sphere
    
    DiMESloc_cart=[DiMESloc_cyl[0]*np.cos(DiMESloc_cyl[1]),
                  DiMESloc_cyl[0]*np.sin(DiMESloc_cyl[1]),DiMESloc_cyl[2]]
    
    xD=np.linspace(DiMESloc_cart[0]-DiMESdims[0],DiMESloc_cart[0]+DiMESdims[0],101)
    yD=np.linspace(DiMESloc_cart[1]-DiMESdims[0],DiMESloc_cart[1]+DiMESdims[0],101)
    
    x1D=np.linspace(DiMESloc_cart[0]-DiMESdims[0],DiMESloc_cart[0]+DiMESdims[0],11)
    y1D=np.linspace(DiMESloc_cart[1]-DiMESdims[0],DiMESloc_cart[1]+DiMESdims[0],11)
    
    xxD,yyD=np.meshgrid(xD,yD,indexing='ij')
    
    arg1=DiMESdims[0]**2-(xxD-DiMESloc_cart[0])**2-(yyD-DiMESloc_cart[1])**2
    arg1[arg1<0]=0
    
    #zsurf=DiMESloc_cart(3)+(DiMESdims(2)/DiMESdims(1))*sqrt(arg1);
    zsurf=DiMESloc_cart[2]-(DiMESdims[0]-DiMESdims[1])+np.sqrt(arg1)
    zsurf[zsurf<DiMESloc_cart[2]]=DiMESloc_cart[2]
    
    rmscale=1.;
    DiMESconflag=np.zeros(np.shape(flagCon));
    DiMESconflag[(xx-DiMESloc_cart[0])**2+(yy-DiMESloc_cart[1])**2+
                (zz-(DiMESloc_cart[2]-(DiMESdims[0]-DiMESdims[1])))**2<(DiMESdims[0]*rmscale)**2]=1;
    
    DiMESdeconfined=np.sum(DiMESconflag,axis=1)

sav=0
if sav==1:
    
    filename='DiMES_impacts_'+run_directory[0]+'.h5'

    NRE_DiMES=int(np.sum(DiMESconflag[-1,:]))

    with h5py.File(filename, "w") as f:
        dset=f.create_dataset('NRE',(1,),dtype='i')
        dset[0]=NRE_DiMES
        dset=f.create_dataset('X',(NRE_DiMES,),dtype='f')
        dset[:]=xx[-1,DiMESconflag[-1,:]>0]
        dset=f.create_dataset('Y',(NRE_DiMES,),dtype='f')
        dset[:]=yy[-1,DiMESconflag[-1,:]>0]
        dset=f.create_dataset('Z',(NRE_DiMES,),dtype='f')
        dset[:]=zz[-1,DiMESconflag[-1,:]>0]
        dset=f.create_dataset('VX',(NRE_DiMES,),dtype='f')
        dset[:]=vx[-1,DiMESconflag[-1,:]>0]
        dset=f.create_dataset('VY',(NRE_DiMES,),dtype='f')
        dset[:]=vy[-1,DiMESconflag[-1,:]>0]
        dset=f.create_dataset('VZ',(NRE_DiMES,),dtype='f')
        dset[:]=vz[-1,DiMESconflag[-1,:]>0]

#%% Save parallel current projection variables

sav=1
if sav==1:
    
    filename='Jpll_projection_vars.h5'

    NRE_DiMES=int(np.sum(DiMESconflag[-1,:]))

    with h5py.File(filename, "w") as f:
        dset=f.create_dataset('NRE',(1,),dtype='i')
        dset[0]=NRE_DiMES
        dset=f.create_dataset('X',(NRE_DiMES,),dtype='f')
        dset[:]=xx[-1,DiMESconflag[-1,:]>0]
        dset=f.create_dataset('Y',(NRE_DiMES,),dtype='f')
        dset[:]=yy[-1,DiMESconflag[-1,:]>0]
        dset=f.create_dataset('Z',(NRE_DiMES,),dtype='f')
        dset[:]=zz[-1,DiMESconflag[-1,:]>0]
        dset=f.create_dataset('VX',(NRE_DiMES,),dtype='f')
        dset[:]=vx[-1,DiMESconflag[-1,:]>0]
        dset=f.create_dataset('VY',(NRE_DiMES,),dtype='f')
        dset[:]=vy[-1,DiMESconflag[-1,:]>0]
        dset=f.create_dataset('VZ',(NRE_DiMES,),dtype='f')
        dset[:]=vz[-1,DiMESconflag[-1,:]>0]

#%% Plotting

#plt.rcParams['text.usetex'] = True
plt.rcParams['contour.negative_linestyle'] = 'solid'

SMALL_SIZE = 12
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=SMALL_SIZE)
plt.rc('figure', titlesize=SMALL_SIZE)

plot_histrm=0
plot_LAC_ParamScaling=0
plot_LAC_Escaling=0
plot_GPUscaling=0
plot_LACbench=0
plot_3Dloc=0
plot_evo=0
plot_orbit=0
plot_histRZ_analytic=0
plot_evoCon=0
plotgrowth=0
plotyorbit=0   
plot_ne=0
plot_Te=0
plot_ephi=0
plot_br=0
plot_bz=0
plot_bphi=0
plot_psip=0
plot_histRZ_ext=0
plot_histKeta=0
plot_histK=0
plot_histeta=0
plot_evoCon_DiMES=0
plot_evo1D=0
plot_fieldm=0
plot_deconloc=0
plot_deconloc1=0
plot_deconloc2=0
plot_deconloc3=0
plot_deconhistphi=0
plot_deconsurf=0
plot_deconhist=0
plot_inc_ang=0
tmpplot=0
plot_psi=0
plot_histRZ_m3dc1=0
plotallangle_fourplot = 0
plot_histtmp = 0
plot_evo0D = 1

timeind_p=0
timeind_g=0

tloss=time[12]
t0=1.594691872596741e+00

need_exp_data=0
if need_exp_data==1:
    filename_ip = '/home/21b/KORC_RUNS/FROM_PERLMUTTER/TEST21/ip177031.txt'

    # Downsampling option: 1 means keep all rows, 2 means take every second row, etc.
    downsample_stride = 1  
    
    # --- Phase 1: Load and Resample ASCII Signal Data ---
    if not os.path.exists(filename_ip):
        raise FileNotFoundError(f"Target data file not found at: {filename_ip}")
    
    # Load rows dynamically using native NumPy streaming
    # Apply stride indexing [::downsample_stride] across rows directly during slice allocation
    raw_data = np.loadtxt(filename_ip)[::downsample_stride]
    
    # Parse columns: Col 1 is time in ms (convert to seconds), Col 2 is Current in Amps
    Ip_time = raw_data[:, 0] / 1e3  # convert ms -> s
    Ip = raw_data[:, 1]             # Current in [A]
    

if need_exp_data==2:
    folder = '/home/21b/KORC_RUNS/FROM_PERLMUTTER/TEST21'
    h5_filename = os.path.join(folder, 'ipprobesf.h5')
    
    # Downsampling option (1 = full data, 2 = every second point, etc.)
    downsample_stride = 10
    
    with h5py.File(h5_filename, 'r') as hf:
        # Define the deep internal path specified by your diagnostic file
        hdf5_path = 'DATA/axes_47589961407024/step_1/args'
            
        target_node = hf[hdf5_path]
        
        # CASE A: If 'args' is a dataset containing a 2D array (e.g., [N x 2])
        if isinstance(target_node, h5py.Dataset):
            raw_data = target_node[::downsample_stride]
            # Adjust indices if your columns are inverted (e.g., Col 0 vs Col 1)
            Ip_time = raw_data[:, 0] / 1e3  # Assuming ms -> seconds conversion
            Ip = raw_data[:, 1]
            
        # CASE B: If 'args' is a group containing separate datasets (e.g., 'time' and 'current')
        elif isinstance(target_node, h5py.Group):
            # Adjust strings inside brackets to match the printed keys if necessary
            Ip_time = target_node['0'][::downsample_stride] / 1e3  
            Ip = target_node['1'][::downsample_stride]
if need_exp_data!=0:
                   # Time base offset scalar

    # Replicating MATLAB's find(timeoffset):tmpfld= bR_inc_2d.copy()
    # We find where target_value > 0. The first match minus 1 maps to index(1) - 1.
    # Using np.where returns an array of indices matching the criteria.
    valid_indices = np.where(Ip_time >= t0)[0]
    
    if len(valid_indices) > 0:
        tindex = valid_indices[0] - 1  # Replicating your index(1) - 1 shift
        
        # Bound protection to prevent an index of -1 if the first element matched
        if tindex < 0:
            tindex = 0
    else:
        # Fallback if the target time is completely out of the bounds of the file
        tindex = len(Ip_time) - 1

if plot_evo0D==1:
    
    tmpfld=Itot.copy()
    tmpfld1=Ipri.copy()
    tmpfld2=Isec.copy()
    
    #Irat=Ip[tindex]/-tmpfld[1]
    Irat=1
    
    fig,ax=plt.subplots()
    
    ax.plot(t0+time[1:]-time[1],-Irat*tmpfld[1:],'-o', label=r'All RE')
    ax.plot(t0+time[1:]-time[1],-Irat*tmpfld1[1:],'-o', label=r'Primary RE')
    ax.plot(t0+time[1:]-time[1],-Irat*tmpfld2[1:],'-o', label=r'Secondary RE')
    #ax.plot(Ip_time[tindex:],Ip[tindex:], label=r'DIIID 177031')
    
    ax.legend(loc='center left', frameon=False, fontsize=12)
    ax.set_xlim([t0,t0+0.010])
    ax.set(xlabel='$t (\\mathrm{s})$')
    ax.grid()
    
    #plt.savefig("RZhist.png", format="png", bbox_inches="tight")
    plt.show()

if plot_histtmp==1:
    
    #tmpfld= vz.copy()
    #tmpfld= bX_inc_2d.copy()
    #tmpfld= bY_inc_2d.copy()
    #tmpfld= bZ_inc_2d.copy()
    #tmpfld= bR_inc_2d.copy()
    #tmpfld= bPHI_inc_2d.copy()
    #tmpfld= np.cos(np.radians(eta_inc_2d.copy()))
    #tmpfld= vx.copy()
    #tmpfld= vy.copy()
    #tmpfld= vmag_inc_2d.copy()  
    #tmpfld=xint.copy()
    #tmpfld= chi_inc_2d.copy()
    #tmpfld=GR[indt_idx,jj_idx]
    tmpfld=zint
    
    mintmp=np.min(tmpfld[(~np.isnan(tmpfld)) & np.isclose(time_inc_2d,tloss)])
    maxtmp=np.max(tmpfld[(~np.isnan(tmpfld)) & np.isclose(time_inc_2d,tloss)])
    
    #mintmp=np.min(tmpfld[(~np.isnan(tmpfld)) & np.isclose(time_inc_1d,tloss)])
    #maxtmp=np.max(tmpfld[(~np.isnan(tmpfld)) & np.isclose(time_inc_1d,tloss)])
    
    incbins=np.linspace(mintmp,maxtmp,14)
    
    H, xedges= np.histogram(tmpfld[(~np.isnan(tmpfld)) & np.isclose(time_inc_2d,tloss)],bins=incbins)\
    #H, xedges= np.histogram(tmpfld[(~np.isnan(tmpfld)) & np.isclose(time_inc_1d,tloss)],bins=incbins)
    
    fig,ax=plt.subplots()
    
    ax.plot(xedges[:-1],H,'-o')
    
    
    ax.set(xlabel='$V_x (\\mathrm{m/s})$', ylabel='$N_{RE}$')
    ax.grid()
    #plt.gca().set_aspect('equal')
    
    #plt.savefig("RZhist.png", format="png", bbox_inches="tight")
    plt.show()

if plot_deconloc3==1:
    
    fig,ax=plt.subplots()
    
    #tmpfld= bX_inc_2d.copy()
    #tmpfld= bY_inc_2d.copy()
    #tmpfld= bZ_inc_2d.copy()
    #tmpfld= bR_inc_2d.copy()
    #tmpfld= bPHI_inc_2d.copy()
    #tmpfld= vx.copy()
    #tmpfld= vy.copy()
    #tmpfld= inc_inc_2d.copy()
    #tmpfld= np.cos(np.radians(eta_inc_2d.copy()))
    #tmpfld= vmag_inc_2d.copy() 
    #tmpfld= chi_inc_2d.copy()
    #tmpfld=xint.copy()
    #tmpfld=GR[indt_idx,jj_idx]
    
    mintmp=np.min(tmpfld[(~np.isnan(tmpfld)) & np.isclose(time_inc_2d,tloss)])
    maxtmp=np.max(tmpfld[(~np.isnan(tmpfld)) & np.isclose(time_inc_2d,tloss)])
    
    #mintmp=np.min(tmpfld[(~np.isnan(tmpfld)) & np.isclose(time_inc_1d,tloss)])
    #maxtmp=np.max(tmpfld[(~np.isnan(tmpfld)) & np.isclose(time_inc_1d,tloss)])
    
    plot_bwr=0
    if plot_bwr==1:
        if maxtmp>np.abs(mintmp):
            mintmp=-maxtmp
        else:
            maxtmp=np.abs(mintmp)
        
    ticklabels=np.linspace(mintmp,maxtmp,7)
    
    if plot_bwr==1:
        sc=ax.scatter(yint[(~np.isnan(tmpfld)) & np.isclose(time_inc_2d,tloss)],zint[(~np.isnan(tmpfld)) & np.isclose(time_inc_2d,tloss)],c=tmpfld[(~np.isnan(tmpfld)) & np.isclose(time_inc_2d,tloss)],alpha=0.8,cmap='bwr',vmin=mintmp,vmax=maxtmp)
    else:
        sc=ax.scatter(yint[(~np.isnan(tmpfld)) & np.isclose(time_inc_2d,tloss)],zint[(~np.isnan(tmpfld)) & np.isclose(time_inc_2d,tloss)],c=tmpfld[(~np.isnan(tmpfld)) & np.isclose(time_inc_2d,tloss)],alpha=0.8,vmin=mintmp,vmax=maxtmp)
 
    #if plot_bwr==1:
    #    sc=ax.scatter(yint_1d[(~np.isnan(tmpfld)) & np.isclose(time_inc_1d,tloss)],zint_1d[(~np.isnan(tmpfld)) & np.isclose(time_inc_1d,tloss)],c=tmpfld[(~np.isnan(tmpfld)) & np.isclose(time_inc_1d,tloss)],alpha=0.8,cmap='bwr',vmin=mintmp,vmax=maxtmp)
    #else:
    #    sc=ax.scatter(yint_1d[(~np.isnan(tmpfld)) & np.isclose(time_inc_1d,tloss)],zint_1d[(~np.isnan(tmpfld)) & np.isclose(time_inc_1d,tloss)],c=tmpfld[(~np.isnan(tmpfld)) & np.isclose(time_inc_1d,tloss)],alpha=0.8,vmin=mintmp,vmax=maxtmp)
    
    
    
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("$V_x({\\rm T})$")


    ax.set(xlabel='$Y\,({\\rm m})$', ylabel='$Z\,({\\rm m})$')
    ax.grid()
    #plt.gca().set_aspect('equal')
    ax.axis([-half_tile,half_tile,-1,1])
    
    plt.show()

if plot_deconloc2==1:
    
    fig,ax=plt.subplots()
    
    plotall=1
    plotpri=0
    plotsec=0
    
    primask_2d  = np.repeat(primary_mask[:, np.newaxis], nsam_chi, axis=1).ravel()
    secmask_2d  = np.repeat(secondary_mask[:, np.newaxis], nsam_chi, axis=1).ravel()
    
    
    if collapse==1:
        if plotall==1:
            #ax.scatter(yint,xint,c='b')
            ax.quiver(yint[(~np.isnan(vx)) & np.isclose(time_inc_2d,tloss)],
                      xint[(~np.isnan(vx)) & np.isclose(time_inc_2d,tloss)],
                      vy[(~np.isnan(vx)) & np.isclose(time_inc_2d,tloss)],
                      vx[(~np.isnan(vx)) & np.isclose(time_inc_2d,tloss)],
                      scale=2e9)
            ax.set(title='All REs')
        if plotpri==1:
            ax.scatter(yint[primask_2d],xint[primask_2d],c='b')
            ax.set(title='Primary REs')
        if plotsec==1:
            ax.scatter(yint[secmask_2d],xint[secmask_2d],c='b')
            ax.set(title='Secondary REs')
    else:
        ax.scatter(PHIint[conlossind==1],Zint[conlossind==1],c='b')

    ax.set(xlabel='$Y({\\rm m})$', ylabel='$X({\\rm m})$')
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.grid()
    #plt.gca().set_aspect('equal')
    #ax.axis([-0.015,0.015,1.015,1.016])
    
    plt.show()

if plotallangle_fourplot == 1:
    
    calchistall = 0
    maxZ = 1.0             # Set to your actual geometric Z-ceiling
    minZ = -1.0            # Set to your actual geometric Z-floor
    max_inc = 1e0          # Limit threshold boundary
    
    incbinedges = np.linspace(0, 1, 50)
    incbinedges_GC = np.linspace(0, 8, 50)
    incbinedges_eta = np.linspace(0, 180, 50)
    incbinedges_chi = np.linspace(0, 360, 50)
    
    
    if calchistall == 1:
        # Replicate your 1D properties into 2D grids (Shape: N_incidents, nsam_chi)
        chi_2d  = chi[indt_idx, jj_idx, :]
        inc_2d  = np.repeat(inc[indt_idx, jj_idx][:, np.newaxis], nsam_chi, axis=1)
        eta_2d  = np.repeat(eta[indt_idx, jj_idx][:, np.newaxis], nsam_chi, axis=1)
        z_2d    = np.repeat(Zint[indt_idx, jj_idx][:, np.newaxis], nsam_chi, axis=1)
        p_2d    = np.repeat(flagPrimary[indt_idx, jj_idx][:, np.newaxis], nsam_chi, axis=1)
        s_2d    = np.repeat(flagSecondary[indt_idx, jj_idx][:, np.newaxis], nsam_chi, axis=1)
        pfc_2d  = np.repeat(PFCinter[indt_idx, jj_idx][:, np.newaxis], nsam_chi, axis=1)
        loss_2d = np.repeat(conlossind[indt_idx, jj_idx][:, np.newaxis], nsam_chi, axis=1)
        
        # Flatten them immediately into 1D sequences to free temporary memory space
        chi_flat  = np.degrees(chi_2d.ravel())
        inc_flat  = inc_2d.ravel()
        eta_flat  = eta_2d.ravel()
        z_flat    = z_2d.ravel()
        p_flat    = p_2d.ravel()
        s_flat    = s_2d.ravel()
        pfc_flat  = pfc_2d.ravel()
        loss_flat = loss_2d.ravel()
        
        # --- 3. STREAMLINED PHYSICS MATH ---
        eta_rad = np.radians(eta_flat)
        inc_clipped = np.clip(inc_flat, -1.0, 1.0)
        
        # Compute incident_angle_all for ONLY your active 1D vector stream
        incident_angle_all = (
            -np.cos(np.arcsin(inc_clipped)) * np.sin(eta_rad) * np.sin(np.radians(chi_flat))
            + inc_flat * np.abs(np.cos(eta_rad))
        )
        
        # ==========================================================================
        # 4. STRICT FILTER MASKS (Forcing ALL invalid slots to NaN)
        # ==========================================================================
        # Define a definitive mask for what constitutes a valid, physical impact event
        # A point is ONLY valid if it stays within Z bounds, is a confirmed limiter strike (pfc),
        # is a tracked transition (loss), and yields a valid geometric forward angle (>0).
        valid_impact_event = (
            (z_flat <= maxZ) & (z_flat >= minZ) & 
            (pfc_flat == 1) & 
            (loss_flat == 1) & 
            (incident_angle_all > 0) & 
            (np.abs(incident_angle_all) <= max_inc)
        )
    
        # Force absolutely everything that fails this condition to NaN
        incident_angle_all[~valid_impact_event] = np.nan
    
        # Generate primary and secondary subsets by ensuring they are valid impacts FIRST
        incident_angle_primary = np.full_like(incident_angle_all, np.nan)
        is_valid_pri = valid_impact_event & (p_flat == 1)
        incident_angle_primary[is_valid_pri] = incident_angle_all[is_valid_pri]
    
        incident_angle_secondary = np.full_like(incident_angle_all, np.nan)
        is_valid_sec = valid_impact_event & (s_flat == 1)
        incident_angle_secondary[is_valid_sec] = incident_angle_all[is_valid_sec]
    
        # ==========================================================================
        # 5. HIGH-SPEED HISTOGRAM COUNTS (Stripping out NaNs completely)
        # ==========================================================================
        # Panel 1: Sin(Theta) counts - Drop all NaNs safely before binning
        hall, _ = np.histogram(incident_angle_all[~np.isnan(incident_angle_all)], bins=incbinedges)
        hpri, _ = np.histogram(incident_angle_primary[~np.isnan(incident_angle_primary)], bins=incbinedges)
        hsec, _ = np.histogram(incident_angle_secondary[~np.isnan(incident_angle_secondary)], bins=incbinedges)
    
        # Panel 2: Guiding center angles (PDF) - Apply the identical masks
        hall_GC, _ = np.histogram(np.degrees(inc_flat[valid_impact_event]), bins=incbinedges_GC, density=True)
        hpri_GC, _ = np.histogram(np.degrees(inc_flat[is_valid_pri]), bins=incbinedges_GC, density=True)
        hsec_GC, _ = np.histogram(np.degrees(inc_flat[is_valid_sec]), bins=incbinedges_GC, density=True)
    
        # Panel 3: Pitch Angles (PDF)
        hall_eta, _ = np.histogram(eta_flat[valid_impact_event], bins=incbinedges_eta, density=True)
        hpri_eta, _ = np.histogram(eta_flat[is_valid_pri], bins=incbinedges_eta, density=True)
        hsec_eta, _ = np.histogram(eta_flat[is_valid_sec], bins=incbinedges_eta, density=True)
    
        # Panel 4: Gyrophases (PDF) - This completely strips out the uninitialized 0.0 values!
        hall_chi, _ = np.histogram(chi_flat[valid_impact_event], bins=incbinedges_chi, density=True)
        hpri_chi, _ = np.histogram(chi_flat[is_valid_pri], bins=incbinedges_chi, density=True)
        hsec_chi, _ = np.histogram(chi_flat[is_valid_sec], bins=incbinedges_chi, density=True)

        
        # Compute bin centers uniformly
        incbins     = (incbinedges[:-1] + incbinedges[1:]) / 2
        incbins_GC  = (incbinedges_GC[:-1] + incbinedges_GC[1:]) / 2
        incbins_eta = (incbinedges_eta[:-1] + incbinedges_eta[1:]) / 2
        incbins_chi = (incbinedges_chi[:-1] + incbinedges_chi[1:]) / 2
        
        # Global normalization scalers
        Nall = np.sum(incident_angle_all > 0)
        dbin = incbinedges[1] - incbinedges[0]

    # --- 3. MATPLOTLIB FIGURE GRAPH GENERATION ---
    # Setup standard 1400x350 proportions (scaled to inches for matplotlib: 14x3.5)
    fig, axs = plt.subplots(1, 4, figsize=(14, 3.5))
    
    # Global typesetting controls to mimic your standard 'Fontsize', 20 setup
    plt.rcParams.update({'font.size': 14, 'axes.linewidth': 1.5})
    
    # ------------------ SUBPLOT 1: SIN(THETA) ------------------
    ax = axs[0]
    norm_factor = Nall * dbin if Nall > 0 else 1.0
    ax.plot(incbins, hall / norm_factor, linewidth=3, color=(0, 0, 1))
    ax.plot(incbins, hpri / norm_factor, linewidth=3, color=(1, 0, 0))
    ax.plot(incbins, hsec / norm_factor, linewidth=3, color=(0, 0.5, 0))
    ax.set_xlim(0, 1)
    ax.set_ylabel(r'$f_{{\rm RE},*}\times N_*/N_{\rm all}$')
    ax.set_xlabel(r'${\rm sin}\,\Theta$')
    
    # Replicate your manual layout positioning for the Panel 1 legend
    ax.set_box_aspect(1) 
    
    # ------------------ SUBPLOT 2: THETA_GC ------------------
    ax = axs[1]
    ax.plot(incbins_GC, hall_GC, linewidth=3, color=(0, 0, 1), label=r'$f_{\rm all}$')
    ax.plot(incbins_GC, hpri_GC, linewidth=3, color=(1, 0, 0), label=r'$f_{\rm init}$')
    ax.plot(incbins_GC, hsec_GC, linewidth=3, color=(0, 0.5, 0), label=r'$f_{\rm sec}$')
    ax.set_xlim(0, 8)
    ax.set_ylabel(r'$f_{{\rm RE},*}$')
    ax.set_xlabel(r'$\theta_{\rm GC}\,(^\circ)$')
    ax.legend(loc='center left', bbox_to_anchor=(0.025, 0.7), frameon=False, fontsize=12)
    ax.set_box_aspect(1)
    
    # ------------------ SUBPLOT 3: ETA PITCH ANGLE ------------------
    ax = axs[2]
    ax.plot(incbins_eta, hall_eta, linewidth=3, color=(0, 0, 1))
    ax.plot(incbins_eta, hpri_eta, linewidth=3, color=(1, 0, 0))
    ax.plot(incbins_eta, hsec_eta, linewidth=3, color=(0, 0.5, 0))
    ax.set_xlim(0, 180)
    ax.set_xticks([0, 45, 90, 135, 180])
    ax.set_xlabel(r'$\eta\,(^\circ)$')
    ax.set_box_aspect(1)
    
    # ------------------ SUBPLOT 4: CHI GYROPHASE ------------------
    ax = axs[3]
    ax.plot(incbins_chi, hall_chi, linewidth=3, color=(0, 0, 1))
    ax.plot(incbins_chi, hpri_chi, linewidth=3, color=(1, 0, 0))
    ax.plot(incbins_chi, hsec_chi, linewidth=3, color=(0, 0.5, 0))
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xlabel(r'$\chi\,(^\circ)$')
    ax.set_box_aspect(1)
            
    # Save the complete asset layout directly to disk matching EPS specifications
    plt.tight_layout()
    plt.show()

if plot_deconhistphi==1:
    
    fig, ax = plt.subplots(figsize=(8, 5))

# Plot a step histogram of the folded incident angles
    ax.hist(
        phi_wrapped, 
        bins=PHIbinedges, 
        edgecolor='black', 
        facecolor='skyblue', 
        alpha=0.7, 
        label='Incident Particles'
    )

    # Label the plot based on the collapse state
    if collapse == 1:
        ax.set_title("Particle Impact Distribution Collapsed Onto One Tile")
        ax.set_xlabel("Folded Poloidal Angle $\phi$ (rad) [0 to $2\pi/n_t$]")
    else:
        ax.set_title("Full Poloidal Impact Distribution Around the Limiter")
        ax.set_xlabel("Poloidal Angle $\phi$ (rad) [0 to $2\pi$]")
    
    ax.set_ylabel("Particle Count")
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend()
    
    plt.tight_layout()
    plt.show()

if plot_deconloc1==1:
    
    fig,ax=plt.subplots()
    
    plotall=0
    plotpri=0
    plotsec=1
    
    if collapse==1:
        if plotall==1:
            ax.scatter(phi_wrapped,Zint[indt_idx,jj_idx],c='b')
            ax.set(title='All REs')
        if plotpri==1:
            ax.scatter(phi_wrapped[primary_mask],Zint[indt_idx[primary_mask],jj_idx[primary_mask]],c='b')
            ax.set(title='Primary REs')
        if plotsec==1:
            ax.scatter(phi_wrapped[secondary_mask],Zint[indt_idx[secondary_mask],jj_idx[secondary_mask]],c='b')
            ax.set(title='Secondary REs')
    else:
        ax.scatter(PHIint[conlossind==1],Zint[conlossind==1],c='b')

    ax.set(xlabel='$\\phi ({\\rm rad})$', ylabel='$Z({\\rm m})$')
    ax.grid()
    #plt.gca().set_aspect('equal')
    ax.axis([-half_tile,half_tile,-1,1])
    
    plt.show()

if plot_histRZ_m3dc1==1:
    
    H, xedges, yedges = np.histogram2d(R[timeind_p,:],zz[timeind_p,:])
    
    fig,ax=plt.subplots()
    
    cmap = plt.cm.viridis
    cmap_modified = cmap.with_extremes(under='white')
    
    ct=ax.pcolormesh(xedges, yedges, H, cmap=cmap_modified,vmin=1)
    
    cb=plt.colorbar(ct)
    cb.set_label('$N_{RE}$', rotation=90)
    
    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    plt.gca().set_aspect('equal')
    
    plt.savefig("RZhist.png", format="png", bbox_inches="tight")
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
        
    nlevs=25
    levs=np.linspace(mintmp,maxtmp,nlevs)
    ticklabels=np.linspace(mintmp,maxtmp,7)
    
    if plot_bwr==1:
        ct=ax.contour(RRm,ZZm,tmpfld,levs,cmap='bwr')
    else:
        ct=ax.contour(RRm,ZZm,tmpfld,levs)
        
    plt.plot(R,zz)
        
    plt.plot(limR,limZ,color='k',linewidth=3)

    cbar=plt.colorbar(ct,ticks=ticklabels,format="%4.1f")
    cbar.set_label('$\psi$', fontsize=12)

    plt.gca().set_aspect('equal')
    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    

    plt.savefig("psi_plot.pdf", format="pdf", bbox_inches="tight")
    plt.show() 

if tmpplot==1:
    fig,ax=plt.subplots()

    ax.plot(np.where(R[timeind_p,flagActive[timeind_p,:]>0]),np.where(R[timeind_p,flagActive[timeind_p,:]>0]),'o')

    ax.grid()
    
    plt.show()

if plot_histrm==1:
    
    rmbin=np.linspace(0,1.2*ar,40)
    
    fig,ax=plt.subplots()

    H,xedges=np.histogram(rm[0,flagActive[0,:]>0],bins=rmbin)
    ax.plot(xedges[:-1],H/xedges[1:],linewidth=2,color='k')
    
    H,xedges=np.histogram(rm[10,flagActive[10,:]>0],bins=rmbin)
    ax.plot(xedges[:-1],H/xedges[1:],linewidth=2,color='r',)
    
    H,xedges=np.histogram(rm[20,flagActive[20,:]>0],bins=rmbin)
    ax.plot(xedges[:-1],H/xedges[1:],linewidth=2,color='b')
    
    H,xedges=np.histogram(rm[30,flagActive[30,:]>0],bins=rmbin)
    ax.plot(xedges[:-1],H/xedges[1:],linewidth=2,color='g')
    
    H,xedges=np.histogram(rm20b3[0,flagActive20b3[0,:]>0],bins=rmbin)
    ax.plot(xedges[:-1],H/xedges[1:],linewidth=2,color='k',linestyle='--')
    
    H,xedges=np.histogram(rm20b3[10,flagActive20b3[10,:]>0],bins=rmbin)
    ax.plot(xedges[:-1],H/xedges[1:],linewidth=2,color='r',linestyle='--')
    
    H,xedges=np.histogram(rm20b3[20,flagActive20b3[20,:]>0],bins=rmbin)
    ax.plot(xedges[:-1],H/xedges[1:],linewidth=2,color='b',linestyle='--')
    
    H,xedges=np.histogram(rm20b3[30,flagActive20b3[30,:]>0],bins=rmbin)
    ax.plot(xedges[:-1],H/xedges[1:],linewidth=2,color='g',linestyle='--')
    
    ax.legend(['$t=0$ No wall','$t=1\,\\mathrm{ms}$','$t=2\,\\mathrm{ms}$',
              '$t=3\,\\mathrm{ms}$','$t=0$ LFS'],loc='lower left',fontsize=13)
    
    #ax.legend([f't= {time[0]:4.1e} s',f't= {time[10]:4.1e} s',f't= {time[20]:4.1e} s',
    #           f't= {time[30]:4.1e} s',f't= {time[40]:4.1e} s'],loc='upper right')
    
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    
    ax.axis([0.6,1.15,-10,400])
    
    ax.set(xlabel='$r_m\,[\\mathrm{m}]$', ylabel='$N_{\\mathrm{RE}}/r_m$')
    ax.grid()
    
    plt.savefig("histrm_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    plt.show()

if plot_LAC_ParamScaling==1:

    fig,(ax1,ax2,ax3)=plt.subplots(1,3,figsize=(12,4))
    
    #Ephi=np.array([2,4,8,16,32])
    #E_CH=0.754
    #E_norm=Ephi/E_CH

    a=[0.33,0.66,1.0,1.33]
    b0=[1.2,2.2,3.2,4.2,8.2,12.2]
    q0=[0.5,1.5,2.5,3.5,4.5]

    tau_ckorc=2.2603950975747852E-003
    tau_c0=0.0016713116113300154

    T20_a_gr =np.array([0.3461,0.2505,0.1998,0.1432])
    T20b_a_gr=np.array([0.2809,0.2131,0.1767,0.1257])
    T20e_a_gr=0.4959

    T20_b_gr =np.array([0.1997,0.1998,0.1931,0.1873,0.2099,0.1954])
    T20b_b_gr=np.array([0.1575,0.1767,0.1944,0.1893,0.1943,0.1996])
    T20e_b_gr=0.4959

    T20_c_gr =np.array([0.1959,0.1998,0.2263,0.2223,0.2161])
    T20b_c_gr=np.array([0.1866,0.1767,0.1570,0.1477,0.1268])
    T20e_c_gr=0.4959


    ax1.plot(a,T20_a_gr*tau_ckorc/tau_c0,'-o',linewidth=2,color='r',markersize=10)
    ax1.plot(a,T20b_a_gr*tau_ckorc/tau_c0,'-o',linewidth=2,color='g',markersize=10)
    #ax1.axhline(y=T20e_a_gr*tau_ckorc/tau_c0, color='k',linestyle=':', linewidth=2)

    ax2.plot(b0,T20_b_gr*tau_ckorc/tau_c0,'-o',linewidth=2,color='r',markersize=10)
    ax2.plot(b0,T20b_b_gr*tau_ckorc/tau_c0,'-o',linewidth=2,color='g',markersize=10)
    #ax2.axhline(y=T20e_b_gr*tau_ckorc/tau_c0, color='k',linestyle=':', linewidth=2)

    ax3.plot(q0,T20_c_gr*tau_ckorc/tau_c0,'-o',linewidth=2,color='r',markersize=10)
    ax3.plot(q0,T20b_c_gr*tau_ckorc/tau_c0,'-o',linewidth=2,color='g',markersize=10)
    #ax3.axhline(y=T20e_c_gr*tau_ckorc/tau_c0, color='k',linestyle=':', linewidth=2)

    #plt.axhline(y=0, color='k', linestyle=':', linewidth=2)

    #ax.set_xscale('log')
    #ax.set_yscale('log')
    
    ax1.set(xlabel='Beam $a\,[\\mathrm{m}]$', ylabel='$\gamma \\tau_c$')
    ax2.set(xlabel='$B_\phi(R_0)\,[\\mathrm{T}]$', ylabel='')
    ax3.set(xlabel='$q_0$', ylabel='')
    
    ax2.tick_params(axis='y', length=0, labelleft=False)
    ax3.tick_params(axis='y', length=0, labelleft=False)
    
    plt.subplots_adjust(wspace=0.05)

    ax1.set_ylim(0.15, 0.5)
    ax2.set_ylim(0.15, 0.5)
    ax3.set_ylim(0.15, 0.5)
    
    ax1.grid()
    ax2.grid()
    ax3.grid()
    
    ax3.legend(['No wall','LFS'],loc='upper right')
    
    plt.savefig("FP_LAC_ParamScaling.png", format="png", bbox_inches="tight")
    #plt.savefig("Orbit_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    plt.show() 

if plot_LAC_Escaling==1:

    fig,ax=plt.subplots(1)
    
    Ephi=np.array([2,4,8,16,32,64,128])
    E_CH=0.754
    E_norm=Ephi/E_CH
    
    tau_ckorc=2.2603950975747852E-003
    tau_c0=0.0016713116113300154

    T20_gr =np.array([0.0198,0.0425,0.1014,0.1998,0.3966,0.5979,0.6734])
    T20a_gr=np.array([0.0149,0.0423,0.0981,0.1988,0.411,0.6102,0.7167])
    T20b_gr=np.array([0.0141,0.0373,0.0748,0.1767,0.3449,0.6165,0.6455])
    T20c_gr=np.array([0.0149,0.0468,0.0946,0.195,0.3856,0.5828,0.6427])
    T20d_gr=np.array([0.0174,0.0433,0.0909,0.1898,0.3957,0.5511,0.7754])
    #T20e_gr=np.array([0.0276,0.0796,0.186,0.3036,0.3597]) #using Te=2000eV
    T20e_gr=np.array([0.04086,0.1080,0.2384,0.4959,0.8705,0.9831,1.110]) #using Te=2eV

    #ax.plot(E_norm,T20e_gr*tau_ckorc/tau_c0,'-o',linewidth=2,color='k',markersize=10)
    #ax.plot(E_norm,T20_gr*tau_ckorc/tau_c0,'-o',linewidth=2,color='r',markersize=10)
    #ax.plot(E_norm,T20a_gr*tau_ckorc/tau_c0,'-o',linewidth=2,color='b',markersize=10)
    #ax.plot(E_norm,T20b_gr*tau_ckorc/tau_c0,'-o',linewidth=2,color='g',markersize=10)
    #ax.plot(E_norm,T20c_gr*tau_ckorc/tau_c0,'-o',linewidth=2,color='c',markersize=10)
    #ax.plot(E_norm,T20d_gr*tau_ckorc/tau_c0,'-o',linewidth=2,color='m',markersize=10)
    
    pmax=5
    
    #ax.plot(E_norm,T20e_gr*T20_gr,'-o',linewidth=2,color='k',markersize=10)
    ax.plot(E_norm[0:pmax],T20_gr[0:pmax]/T20_gr[0:pmax],'-o',linewidth=2,color='r',markersize=10)
    ax.plot(E_norm[0:pmax],T20a_gr[0:pmax]/T20_gr[0:pmax],'-o',linewidth=2,color='b',markersize=10)
    ax.plot(E_norm[0:pmax],T20b_gr[0:pmax]/T20_gr[0:pmax],'-o',linewidth=2,color='g',markersize=10)
    ax.plot(E_norm[0:pmax],T20c_gr[0:pmax]/T20_gr[0:pmax],'-o',linewidth=2,color='c',markersize=10)
    ax.plot(E_norm[0:pmax],T20d_gr[0:pmax]/T20_gr[0:pmax],'-o',linewidth=2,color='m',markersize=10)
    
    
    #plt.axhline(y=0, color='k', linestyle=':', linewidth=2)

    ax.set_xscale('log')
    #ax.set_yscale('log')
    
    #ax.set(xlabel='$E/E_{\\mathrm{CH}}$', ylabel='$\gamma \\tau_c$')
    
    ax.set(xlabel='$E/E_{\\mathrm{CH}}$', ylabel='$\gamma/\gamma_{\\rm no\,wall}$')
    
    #ax.axis([0.e-6,1e-5,10.7,11.2])
    
    ax.grid()
    
    #ax.legend(['FP','No wall','HFS','LFS','Top','Bottom'],loc='upper left',ncol=2)
    ax.legend(['No wall','HFS','LFS','Top','Bottom'],loc='lower right',ncol=2)
    
    plt.savefig("FP_LAC_Escaling.png", format="png", bbox_inches="tight")
    #plt.savefig("Orbit_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    plt.show() 

if plot_GPUscaling==1:

    fig,ax=plt.subplots(1)
    
    nRE=np.array([10**2,10**3,10**4,10**5,10**6,10**7,10**8])
    
    orbit_time=np.array([0.3477870,0.3611510,0.3603720,0.3895470,0.6593030,3.370124,30.24167])
    orbit_time0=0.2199159472500583
    
    coll_time=np.array([3.956,6.539,37.40,348.7,348.7*10,348.7*100,348.7*1000])
    
    LAC_time=np.array([14.65,30.07,115.4,943.8,943.8*10,943.8*100,943.8*1000])
    
    orbit_time=orbit_time
    coll_time=coll_time
    LAC_time=LAC_time
    
    gpuREs=108*64*32

    ax.plot(nRE,orbit_time,'-o',linewidth=2,color='r',markersize=10)
    ax.plot(nRE,coll_time,'-o',linewidth=2,color='b',markersize=10)
    ax.plot(nRE,LAC_time,'-o',linewidth=2,color='g',markersize=10)

    plt.axvline(x=gpuREs, color='k', linestyle=':', linewidth=2)

    ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.set(xlabel='$N_{\\mathrm{RE}}$', ylabel='Runtime $(\\mathrm{s})$')
    
    #ax.axis([0.5e2,2.e8,0.75,110])
    
    ax.grid()
    
    ax.legend(['Orbit','Coulomb','Large-angle'],loc='upper left')
    
    plt.savefig("GPU_scaling.png", format="png", bbox_inches="tight")
    #plt.savefig("Orbit_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    plt.show() 

if plot_LACbench==1:

    fig,ax=plt.subplots(1)
    
    ne=1e20
    Clog=20

    tau_c=4*np.pi*ep0**2*me**2*c**3/(ne*qe**4*Clog)
    
    E_CH=me*c/(qe*tau_c)

    E_norm=np.array([2.5,2.75,3,3.5,4,5])
    gr5=np.array([-4.27,-0.964,0.645,2.20,3.08,4.59]) #from FP_LAC on Cori
    
    grGPC=np.array([-0.073,-0.021,0.0081,0.035,0.049,0.074])

    ax.plot(E_norm,gr5*tau_c,'-o',linewidth=2,color='r',markersize=10)
    ax.plot(E_norm,grGPC,'-*',linewidth=2,color='b',markersize=10)
    
    plt.axhline(y=0, color='k', linestyle=':', linewidth=2)

    #ax.set_xscale('log')
    #ax.set_yscale('log')
    
    ax.set(xlabel='$E/E_{\\mathrm{CH}}$', ylabel='$\gamma \\tau_c$')
    
    #ax.axis([0.e-6,1e-5,10.7,11.2])
    
    ax.grid()
    
    ax.legend(['Cori CPU','Perlmutter GPU'],loc='lower right')
    
    plt.savefig("FP_LAC_benchmark.png", format="png", bbox_inches="tight")
    #plt.savefig("Orbit_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    plt.show()   

if plot_3Dloc==1:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')   
    ax.view_init(elev=90, azim=0, roll=0)

    nRE=1000
    
    ax.scatter(xx[timeind_p,:nRE],yy[timeind_p,:nRE],zz[timeind_p,:nRE],c='b',marker='.')

    ax.set(xlabel='$X (m)$',ylabel='$Y (m)$',zlabel='$Z (m)$')
    plt.gca().set_aspect('equal')
    
    plt.savefig("3Dloc", format="png", bbox_inches="tight")
    plt.show()

if plot_evo==1:

    fig,ax=plt.subplots(1)
    
    tmp_time=time.copy()
    tmp_R=R.copy()
    tmp_PHI=PHI.copy()
    tmp_Z=zz.copy()
    #tmp=eta.copy()
    tmp=g.copy()

    #ax.plot(tmp_R,tmp_Z,'o-')
    ax.plot(tmp_time,tmp,'.-')

    #ax.set_xscale('log')
    #ax.set_yscale('log')
    
    ax.set(xlabel='$t (\\mathrm{s})$', ylabel='$Z (\\mathrm{m})$')
    
    #ax.axis([0.e-6,1e-5,10.7,11.2])
    
    ax.grid()
    
    plt.savefig("Evo_"+"TEST8"+".png", format="png", bbox_inches="tight")
    #plt.savefig("Orbit_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    plt.show()   

if plot_orbit==1:

    fig,ax=plt.subplots(2,2)
    
    tmp_time=time.copy()
    
    singleorbit=1
    if singleorbit==1:
        pind=0
        
        tmp_R=R[:,pind].copy()
        tmp_PHI=PHI[:,pind].copy()
        tmp_Z=zz[:,pind].copy()
    else:   
        tmp_R=R.copy()
        tmp_PHI=PHI.copy()
        tmp_Z=zz.copy()

    #ax.plot(tmp_R,tmp_Z,'o-')
    ax[0,0].plot(tmp_time,tmp_Z,'.-')
    ax[0,1].plot(tmp_time,tmp_R,'.-')
    ax[1,0].plot(tmp_time,tmp_PHI,'.-')
    ax[1,1].plot(tmp_R,tmp_Z,'.-')
    if singleorbit==1:
        ax[1,1].plot(tmp_R[0],tmp_Z[0],'ro')

    #ax.set_xscale('log')
    #ax.set_yscale('log')
    
    ax[0,0].set(xlabel='$t (\\mathrm{s})$', ylabel='$Z (\\mathrm{m})$')
    ax[0,1].set(xlabel='$t (\\mathrm{s})$', ylabel='$R (\\mathrm{m})$')
    ax[1,0].set(xlabel='$t (\\mathrm{s})$', ylabel='$\phi (\\mathrm{rad})$')
    ax[1,1].set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    
    ax[0,0].grid()
    ax[0,1].grid()
    ax[1,0].grid()
    ax[1,1].grid()
    
    ax[1,1].set_aspect('equal')
    
    plt.subplots_adjust(wspace=0.6, hspace=0.4)
    
    plt.savefig("Orbit_"+"TEST8"+".png", format="png", bbox_inches="tight")
    #plt.savefig("Orbit_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    plt.show()   

if plot_histRZ_analytic==1:

    tmp=np.where

    #H, xedges, yedges = np.histogram2d(R[timeind_p,flagCon[timeind_p,:]<1], 
    #                   zz[timeind_p,flagCon[timeind_p,:]<1],bins=20)

    #H, xedges, yedges = np.histogram2d(R[timeind_p,flagCol[timeind_p,:]>0], 
    #                   zz[timeind_p,flagCol[timeind_p,:]>0],bins=20)

    H, xedges, yedges = np.histogram2d(R[timeind_p,flagActive[timeind_p,:]>0], 
                      zz[timeind_p,flagActive[timeind_p,:]>0],bins=20)

    
    fig,ax=plt.subplots()
    
    cmap = plt.cm.viridis

    # Create a modified colormap with 'white' as the "under" color
    cmap_modified = cmap.with_extremes(under='white')
        

    ct=ax.pcolormesh(xedges, yedges, H.T, cmap=cmap_modified,vmin=1)

    ax.axis([0,5,-3,3])
    
    cb=plt.colorbar(ct)
    

    cb.set_label('$N_{RE}$', rotation=90)

    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    plt.gca().set_aspect('equal')
    
    ax.set(title=f't= {time[timeind_p]:4.1e} s')
    
    plt.savefig("RZhist"+filename[0]+".png", format="png", bbox_inches="tight")
    plt.show()
    

if plot_evoCon==1:

    fig,ax=plt.subplots()
    
    
    tmp_time=time.copy()
    tmp_tot=Total.copy()
    tmp_act=Active.copy()
    tmp_confined=confined.copy()
    tmp_Energetic=Energetic.copy()
    
    ax.plot(tmp_time,tmp_tot,'-')
    ax.plot(tmp_time,tmp_confined,'-')
    ax.plot(tmp_time,tmp_Energetic,'-')

    ax.legend(['Total','Confined','Energetic'])

    #ax.set_xscale('log')
    #ax.set_yscale('log')
    
    ax.set(xlabel='$t\,({\\rm s})$', ylabel='$N_{RE}$')
    ax.grid()
    
    #plt.savefig("NREevo_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    plt.show()    

if plotgrowth==1:
    
    Clog0=20
    ne=1e21
    
    gammac=ne*abs(qe**4)*Clog0/(4*np.pi*ep0**2)
    #%collisional gamma
    rel_col_freq=gammac/(me**2*c**3)
    #%collision frequency based off of relativistic velocity
    E_CH=rel_col_freq*me*c/abs(qe)
    tau_c0=me**2*c**3/gammac
    rad_damp_time=6*np.pi*ep0*(me*c)**3/(qe**4*B0**2)
    
    tmp_time=time.copy()
    tmp_NREtot=Total.copy()
    tmp_NREact=Active.copy()
    
    tmp_time=np.trim_zeros(tmp_time,trim='b')
    tmp_NREtot=np.trim_zeros(tmp_NREtot,trim='b')
    tmp_NREact=np.trim_zeros(tmp_NREact,trim='b')
    
    iffit=1
    
    if iffit==1:
        indrange0=int(np.round(np.shape(tmp_time)[0]*.75))
        
        def exp_func(x,b,c):
            return b*np.exp(c*x)
        
        popt_tot, pcov = sci.optimize.curve_fit(exp_func, tmp_time[indrange0:], tmp_NREtot[indrange0:], p0=[tmp_NREtot[0],1])
        popt_act, pcov = sci.optimize.curve_fit(exp_func, tmp_time[indrange0:], tmp_NREact[indrange0:], p0=[tmp_NREact[0],1])
        
        b_fit_tot,c_fit_tot=popt_tot
        b_fit_act,c_fit_act=popt_act

        print(b_fit_act,c_fit_act*tau_c0)
    
    fig,ax=plt.subplots()
    
    ax.plot(tmp_time,tmp_NREtot,'o-')
    ax.plot(tmp_time,tmp_NREact,'o-')
    
    if iffit==1:
        ax.plot(tmp_time[indrange0:],b_fit_tot*np.exp(c_fit_tot*tmp_time[indrange0:]))
        ax.plot(tmp_time[indrange0:],b_fit_act*np.exp(c_fit_act*tmp_time[indrange0:]))

    #ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.text(0.05,.85,r"$E/E_{\rm CH}=$"+str(np.round(E0/E_CH,2))+"\n"+r"$\gamma \tau_c=$"+str(np.round(c_fit_act*tau_c0,4)),transform=ax.transAxes)
    
    ax.set(xlabel='$t\,({\\rm s})$', ylabel='$N_{RE}$')
    ax.grid()
    
    plt.savefig("NREevo_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    #plt.show()

if plot_inc_ang==1:

    vz_inc=1
    bz_inc=0    

    if bz_inc==1:
        zinc_ang=np.rad2deg(np.arcsin(-bZ/bmag))    
    elif vz_inc==1:
        zinc_ang=np.rad2deg(np.arcsin(-vz/vmag))
    
    etabin=np.linspace(0,10,num=50)
    
    fig,ax=plt.subplots()
    
    plotden=0
    
    if plotden==0:
        H,xedges=np.histogram(zinc_ang[timeind_p,flagCon[timeind_p]<1],bins=etabin)
    else:
        H,xedges=np.histogram(zinc_ang[timeind_p,flagCon[timeind_p]<1],bins=etabin,density=True)

    ax.plot(xedges[:-1],H)
    
    
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    
    ax.set(xlabel='$\\theta_{\\rm inc} (^{\circ})$', ylabel='$N_{RE}$')
    ax.grid()
    
    plt.savefig("zinc_ang_hist_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    plt.show()

if plot_deconhist==1:

    theta=np.linspace(0,2*np.pi,360)
    RDiMES=np.sqrt(DiMESdims[0]**2-(DiMESdims[0]-DiMESdims[1])**2)  
    
    xDiMES=DiMESloc_cart[0]+RDiMES*np.cos(theta)
    yDiMES=DiMESloc_cart[1]+RDiMES*np.sin(theta)
    
    xbin=np.linspace(min(xD),max(xD),num=20)
    ybin=np.linspace(min(yD),max(yD),num=20)
    
    xxbin,yybin=np.meshgrid(xbin,ybin,indexing='ij')
    
    plotden=0
    
    if plotden==0:
        H, xedges, yedges = np.histogram2d(xx[timeind_p,DiMESconflag[timeind_p]>0], 
                                      yy[timeind_p,DiMESconflag[timeind_p]>0], bins=(xbin, ybin))
    else:

        H, xedges, yedges = np.histogram2d(xx[timeind_p,DiMESconflag[timeind_p]>0], 
                                      yy[timeind_p,DiMESconflag[timeind_p]>0], bins=(xbin, ybin),density=True)

    fig,ax=plt.subplots()
    
    cmap = plt.cm.viridis

    # Create a modified colormap with 'white' as the "under" color
    cmap_modified = cmap.with_extremes(under='white')
        

    if plotden==1:
        ct=ax.pcolormesh(xxbin, yybin, H, cmap=cmap_modified,vmin=0.99*np.min(H[H>0]))
    else:
        ct=ax.pcolormesh(xxbin, yybin, H, cmap=cmap_modified,vmin=1)
    
    ax.plot(xDiMES,yDiMES,'k-')
    
    #plt.annotate('$1$',xytext=(-1.26,0.76),xy=(-1.26+0.015*np.cos(DiMESloc_cyl[1]),0.76+0.015*np.sin(DiMESloc_cyl[1])),        
    #        arrowprops=dict(
    #        facecolor='red',  # Color of the arrow
    #        width=2,          # Width of the arrow shaft
    #        headwidth=10,     # Width of the arrow head
    #        headlength=15     # Length of the arrow head
    #    ))
    
    #plt.annotate('$1$',xytext=(-1.26,0.76),xy=(-1.26-0.015*np.sin(DiMESloc_cyl[1]),0.76+0.015*np.cos(DiMESloc_cyl[1])),        
    #        arrowprops=dict(
    #        facecolor='blue',  # Color of the arrow
    #        width=2,          # Width of the arrow shaft
    #        headwidth=10,     # Width of the arrow head
    #        headlength=15     # Length of the arrow head
    #    ))

    import matplotlib.patches as mpatches    
    arr = mpatches.FancyArrowPatch((-1.255,0.765), (-1.255+0.015*np.cos(DiMESloc_cyl[1]),0.765+0.015*np.sin(DiMESloc_cyl[1])),
                              arrowstyle='->,head_width=.15', mutation_scale=20)
    ax.add_patch(arr)
    ax.annotate("$\\hat{R}$", (.65, .45), xycoords=arr, ha='center', va='bottom')
    
    arr = mpatches.FancyArrowPatch((-1.255,0.765), (-1.255-0.015*np.sin(DiMESloc_cyl[1]),0.765+0.015*np.cos(DiMESloc_cyl[1])),
                              arrowstyle='->,head_width=.15', mutation_scale=20)
    ax.add_patch(arr)
    ax.annotate("$\\hat{\phi}$", (.85, .75), xycoords=arr, ha='center', va='top')
    
    cb=plt.colorbar(ct)
    
    if plotden==0:
        cb.set_label('$N_{RE}$', rotation=90)
    else:
        cb.set_label('$f_{RE}$', rotation=90)
    
    ax.set(xlabel='$X (\\mathrm{m})$', ylabel='$Y (\\mathrm{m})$')
    ax.grid()
    plt.gca().set_aspect('equal')
    
    
    plt.savefig("decon_hist_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    plt.show()    

if plot_deconsurf==1:
    
    Plim=np.linspace(0,2*np.pi,360)
    
    RRlim,PPlim=np.meshgrid(Rlim,Plim,indexing='ij')
    
    XXlim=RRlim*np.cos(PPlim)
    YYlim=RRlim*np.sin(PPlim)
    ZZlim=np.zeros(np.shape(XXlim))
    for ii in range(0,np.shape(XXlim)[1]):
        ZZlim[:,ii]=Zlim[:,0]
        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')   
    ax.view_init(elev=0, azim=-90, roll=0)
    
    #vacuum vessel plotting options
    clamshell=0
    noinner=1
    
    DiMESfocus=1
    
    #if DiMESfocus==0:
    #    if clamshell==1:
    #        ax.plot_surface(XXlim[44:,0:180], YYlim[44:,0:180], ZZlim[44:,0:180],color='grey')
    #    elif noinner==1:
    #        ax.plot_surface(XXlim[44:75,0:180], YYlim[44:75,0:180], ZZlim[44:75,0:180],color='grey')    
    #    else:
    #        ax.plot_surface(XXlim[:,0:180], YYlim[:,0:180], ZZlim[:,0:180],color='grey')
    #else:
    #ax.plot_surface(XXlim[44:75,140:160], YYlim[44:75,140:160], ZZlim[44:75,140:160],color='grey',alpha=0.5)
        
    ax.plot_surface(xxD, yyD, zsurf,color='g',zorder=1)
    
    DiMES_decon=1
    
    if DiMES_decon==1:
        #ax.scatter(xx[timeind_p,flagCon[timeind_p]<1],yy[timeind_p,flagCon[timeind_p]<1],
        #           zz[timeind_p,flagCon[timeind_p]<1],c='b')
        ax.scatter(xx[timeind_p,DiMESconflag[timeind_p]>0],yy[timeind_p,DiMESconflag[timeind_p]>0],
                  zz[timeind_p,DiMESconflag[timeind_p]>0],color='r',zorder=4)

    else:
        #ax.scatter(xx[timeind_p,flagCon[timeind_p]<1],yy[timeind_p,flagCon[timeind_p]<1],
        #           zz[timeind_p,flagCon[timeind_p]<1],c='b')
        ax.scatter(xx[timeind_p,:],yy[timeind_p,:],zz[timeind_p,:],c='b')
    
    
    
    if DiMESfocus==1:
        ax.axis([np.min(xD),np.max(xD),np.min(yD),np.max(yD),DiMESloc_cyl[2],DiMESloc_cyl[2]+DiMESdims[1]])
    
    ax.set(xlabel='$X (m)$',ylabel='$Y (m)$',zlabel='$Z (m)$')
    plt.gca().set_aspect('equal')
    
    plt.savefig("decon_surface_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    plt.show()
    

if plot_deconloc==1:
    
    fig,ax=plt.subplots()
    
    DiMES_decon=1
    
    if DiMES_decon==1:
        ax.scatter(PHI[timeind_p,flagCon[timeind_p]<1],theta_geo[timeind_p,flagCon[timeind_p]<1],c='b')
        ax.scatter(PHI[timeind_p,DiMESconflag[timeind_p]>0],theta_geo[timeind_p,DiMESconflag[timeind_p]>0],c='r')

    else:
        ax.scatter(PHI[timeind_p,flagCon[timeind_p]<1],theta_geo[timeind_p,flagCon[timeind_p]<1],c='b')

    ax.set(xlabel='$\\phi ({\\rm rad})$', ylabel='$\\theta_{\\rm geom} ({\\rm rad})$')
    ax.grid()
    #plt.gca().set_aspect('equal')
    #ax.axis([-np.pi,np.pi,-np.pi,np.pi])
    ax.axis([-np.pi,np.pi,-1.8,-1.2])
    
    plt.savefig("decon_loc_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    plt.show()


if plot_histKeta==1:
    
    Kbin10=np.linspace(np.log10(np.min(K[~np.isnan(K)])),np.log10(np.max(K[~np.isnan(K)])),num=50)
    Kbin=10**Kbin10
    #Kbin=np.linspace(np.min(K),np.max(K),num=50)
    
    etabin=np.linspace(0,180,num=50)
    
    KKbin,eetabin=np.meshgrid(Kbin,etabin,indexing='ij')
    
    plotden=0
    plotconfined=1
    plotdeconfined=0
    
    if plotden==0:
        if plotconfined==1:
            H, xedges, yedges = np.histogram2d(K[timeind_p,flagCon[timeind_p]>0], 
                                              eta[timeind_p,flagCon[timeind_p]>0], bins=(Kbin, etabin))
        elif plotdeconfined==1:
            H, xedges, yedges = np.histogram2d(K[timeind_p,flagCon[timeind_p]<1], 
                                              eta[timeind_p,flagCon[timeind_p]<1], bins=(Kbin, etabin))
        else:
            H, xedges, yedges = np.histogram2d(K[timeind_p], eta[timeind_p], bins=(Kbin, etabin))
            
    else:
        if plotconfined==1:
            H, xedges, yedges = np.histogram2d(K[timeind_p,flagCon[timeind_p]>0], 
                                              eta[timeind_p,flagCon[timeind_p]>0], bins=(Kbin, etabin),density=True)
        else:
            H, xedges, yedges = np.histogram2d(K[timeind_p], eta[timeind_p], bins=(Kbin, etabin),density=True)
    
    fig,ax=plt.subplots() 
    
    cmap = plt.cm.viridis

    # Create a modified colormap with 'white' as the "under" color
    cmap_modified = cmap.with_extremes(under='white')
        
    if plotden==1:
        ct=ax.pcolormesh(KKbin, eetabin, H, cmap=cmap_modified,vmin=0.99*np.min(H[H>0]))
    else:
        ct=ax.pcolormesh(KKbin, eetabin, H, cmap=cmap_modified,vmin=1)
    #ax.contour(Rg,Zg,-PSIp,25,colors='black',linewidths=0.5)
    #ax.contour(Rg,Zg,FLAG,[.5],colors='black')

    
    cb=plt.colorbar(ct)
    if plotden==0:
        cb.set_label('$N_{RE}$', rotation=90)
    else:
        cb.set_label('$f_{RE}$', rotation=90)
    
    ax.set(title=f't= {time[timeind_p]:4.1e} s')
    
    ax.set(xlabel='$\\mathcal{K} (\\mathrm{eV})$', ylabel='$\eta (^\circ)$')
    ax.grid()
    #plt.gca().set_aspect('equal')
    ax.set_xscale('log')
    
    plt.savefig("Ketahist_DIIID_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    plt.show()
    
if plot_histK==1:
    
    plotallpri=1
    plotactpri=0
    plotdecpri=0
    plotthepri=0 #no thermalized primaries at t=1.6007s
    
    plotallsec=0
    plotactsec=0
    plotdecsec=0
    plotthesec=0
    
    if plotallpri==1:
        Kbin10=np.linspace(np.log10(np.min(K[flagPrimary>0])),np.log10(np.max(K[flagPrimary>0])),num=50)
    if plotactpri==1:
        Kbin10=np.linspace(np.log10(np.min(K[(flagPrimary>0) & (flagActive>0)])),np.log10(np.max(K[(flagPrimary>0) & (flagActive>0)])),num=50)
    if plotdecpri==1:
        Kbin10=np.linspace(np.log10(np.min(K[(flagPrimary>0) & (flagDecon>0)])),np.log10(np.max(K[(flagPrimary>0) & (flagDecon>0)])),num=50)
    if plotthepri==1:
        Kbin10=np.linspace(np.log10(np.min(K[(flagPrimary>0) & (flagTherm>0)])),np.log10(np.max(K[(flagPrimary>0) & (flagTherm>0)])),num=50)
        
    if plotallsec==1:
        Kbin10=np.linspace(np.log10(np.min(K[flagSecondary>0])),np.log10(np.max(K[flagSecondary>0])),num=50)
    if plotactsec==1:
        Kbin10=np.linspace(np.log10(np.min(K[(flagSecondary>0) & (flagActive>0)])),np.log10(np.max(K[(flagSecondary>0) & (flagActive>0)])),num=50)
    if plotdecsec==1:
        Kbin10=np.linspace(np.log10(np.min(K[(flagSecondary>0) & (flagDecon>0) & ~(flagTherm>0)])),np.log10(np.max(K[(flagSecondary>0) & (flagDecon>0) & ~(flagTherm>0)])),num=50)
    if plotthesec==1:
        Kbin10=np.linspace(np.log10(np.min(K[(flagSecondary>0) & (flagTherm>0) & ~(flagDecon>0)])),np.log10(np.max(K[(flagSecondary>0) & (flagTherm>0) & ~(flagDecon>0)])),num=50)
        
    Kbin=10**Kbin10

    fig,ax=plt.subplots()
    
    legends=[f't= {1.594691872596741e+00+time[ii]:4.4e} s' for ii in range(0,num_snapshots-1)]

    for ii in range(1,np.shape(K)[0]):
        if plotallpri==1:
            H,xedges=np.histogram(K[ii,flagPrimary[ii,:]>0],bins=Kbin)
        if plotactpri==1:
            H,xedges=np.histogram(K[ii,(flagPrimary[ii,:]>0) & (flagActive[ii,:]>0)],bins=Kbin)
        if plotdecpri==1:
            H,xedges=np.histogram(K[ii,(flagPrimary[ii,:]>0) & (flagDecon[ii,:]>0)],bins=Kbin)
        if plotthepri==1:
            H,xedges=np.histogram(K[ii,(flagPrimary[ii,:]>0) & (flagTherm[ii,:]>0)],bins=Kbin)
            
        if plotallsec==1:
            H,xedges=np.histogram(K[ii,flagSecondary[ii,:]>0],bins=Kbin)
        if plotactsec==1:
            H,xedges=np.histogram(K[ii,(flagSecondary[ii,:]>0) & (flagActive[ii,:]>0)],bins=Kbin)
        if plotdecsec==1:
            H,xedges=np.histogram(K[ii,(flagSecondary[ii,:]>0) & (flagDecon[ii,:]>0) & ~(flagTherm[ii,:]>0)],bins=Kbin)
        if plotthesec==1:
            H,xedges=np.histogram(K[ii,(flagSecondary[ii,:]>0) & (flagTherm[ii,:]>0) & ~(flagDecon[ii,:]>0)],bins=Kbin)
    
        ax.plot(xedges[:-1],H,'-o')
    
        #ax.plot(EfitVc*10**6,fM[3,:]/fMnorm[3])
        
        #ax.legend(['Sampled','Target'])
        
        if ii==1:
          ax.set_xscale('log')
          #ax.set_yscale('log')
        
          ax.set(xlabel='$\\mathcal{K} (\\mathrm{eV})$', ylabel='$N_{\\mathrm{RE}}$')

          ax.grid()
        
        #plt.savefig("Khist_DIIID.png", format="png", bbox_inches="tight")
    
    ax.legend(legends,ncol=2,fontsize=8)
    
    if plotallpri==1:
      ax.set(title='All Primary REs')
    if plotactpri==1:
       ax.set(title='Active Primary REs')
    if plotdecpri==1:
        ax.set(title='Deconfined Primary REs')
    if plotthepri==1:
        ax.set(title='Thermalized Primary REs')
        
    if plotallsec==1:
      ax.set(title='All Secondary REs')
    if plotactsec==1:
       ax.set(title='Active Secondary REs')
    if plotdecsec==1:
        ax.set(title='Deconfined Secondary REs')
    if plotthesec==1:
        ax.set(title='Thermalized Secondary REs')
    
    plt.show()


    
if plot_histeta==1:
    
    plotall=0
    plotact=0
    plotdec=1
    plotthe=0
    
    plotallpri=0
    plotactpri=0
    plotdecpri=0
    plotthepri=0 #no thermalized primaries at t=1.6007s
    
    plotallsec=0
    plotactsec=0
    plotdecsec=0
    plotthesec=0
    
    
    tmpfld=eta.copy()
    
    tmpfld=np.where(tmpfld>90,180-tmpfld,tmpfld)
    
    etabin=np.linspace(0,90,50)

    fig,ax=plt.subplots()
    
    legends=[f't= {1.594691872596741e+00+time[ii]:4.4e} s' for ii in range(0,num_snapshots-1)]

    for ii in range(1,np.shape(K)[0]):
        if plotall==1:
            H,xedges=np.histogram(tmpfld[ii],bins=etabin)
        if plotact==1:
            H,xedges=np.histogram(tmpfld[ii, (flagActive[ii,:]>0)],bins=etabin)
        if plotdec==1:
            H,xedges=np.histogram(tmpfld[ii, (flagDecon[ii,:]>0)],bins=etabin)
        if plotthe==1:
            H,xedges=np.histogram(tmpfld[ii, (flagTherm[ii,:]>0)],bins=etabin)
        
        if plotallpri==1:
            H,xedges=np.histogram(tmpfld[ii,flagPrimary[ii,:]>0],bins=etabin)
        if plotactpri==1:
            H,xedges=np.histogram(tmpfld[ii,(flagPrimary[ii,:]>0) & (flagActive[ii,:]>0)],bins=etabin)
        if plotdecpri==1:
            H,xedges=np.histogram(tmpfld[ii,(flagPrimary[ii,:]>0) & (flagDecon[ii,:]>0)],bins=etabin)
        if plotthepri==1:
            H,xedges=np.histogram(tmpfld[ii,(flagPrimary[ii,:]>0) & (flagTherm[ii,:]>0)],bins=etabin)
            
        if plotallsec==1:
            H,xedges=np.histogram(tmpfld[ii,flagSecondary[ii,:]>0],bins=etabin)
        if plotactsec==1:
            H,xedges=np.histogram(tmpfld[ii,(flagSecondary[ii,:]>0) & (flagActive[ii,:]>0)],bins=etabin)
        if plotdecsec==1:
            H,xedges=np.histogram(tmpfld[ii,(flagSecondary[ii,:]>0) & (flagDecon[ii,:]>0) & ~(flagTherm[ii,:]>0)],bins=etabin)
        if plotthesec==1:
            H,xedges=np.histogram(tmpfld[ii,(flagSecondary[ii,:]>0) & (flagTherm[ii,:]>0) & ~(flagDecon[ii,:]>0)],bins=etabin)
    
        ax.plot(xedges[:-1],H,'-o')
    
        #ax.plot(EfitVc*10**6,fM[3,:]/fMnorm[3])
        
        #ax.legend(['Sampled','Target'])
        
        if ii==1:
          #ax.set_xscale('log')
          #ax.set_yscale('log')
        
          ax.set(xlabel='$\\eta (^\\circ)$', ylabel='$N_{\\mathrm{RE}}$')
          ax.xaxis.set_major_locator(MaxNLocator(nbins=7))

          ax.grid()
        
        #plt.savefig("Khist_DIIID.png", format="png", bbox_inches="tight")
    
    ax.legend(legends,ncol=2,fontsize=8)
    
    if plotallpri==1:
      ax.set(title='All Primary REs')
    if plotactpri==1:
       ax.set(title='Active Primary REs')
    if plotdecpri==1:
        ax.set(title='Deconfined Primary REs')
    if plotthepri==1:
        ax.set(title='Thermalized Primary REs')
        
    if plotallsec==1:
      ax.set(title='All Secondary REs')
    if plotactsec==1:
       ax.set(title='Active Secondary REs')
    if plotdecsec==1:
        ax.set(title='Deconfined Secondary REs')
    if plotthesec==1:
        ax.set(title='Thermalized Secondary REs')
    
    plt.show()

if plot_evo1D==1:
    fig,ax=plt.subplots()
    
    tmp=bR
    
    ax.plot(time, tmp)
    
    ax.set(xlabel='$t (\\mathrm{s})$', ylabel='$\phi (\\mathrm{rad})$')
    ax.grid()
    
    plt.savefig("EvoCon_DIIID_AORSA.pdf", format="pdf", bbox_inches="tight")
    plt.show()   

if plot_fieldm==1:
    fig,ax=plt.subplots()
    
    tmpfld=PSIPm.copy()
    
    #nemin=np.min(-tmpfld)
    #nemax=np.max(-tmpfld)
    
    
    ct=ax.contourf(Rm,Zm,tmpfld)
    #ax.contour(Rg,Zg,-PSIp[:,timeind_g,:],25,colors='black',linewidths=0.5)
    #ax.contour(Rg,Zg,FLAG[:,timeind_g,:],[.5],colors='black')
    
    #sc=ax.scatter(R,zz,c=-psiP,s=5,cmap=ct.cmap,vmin=nemin,vmax=nemax,edgecolor='black',linewidth=0.25)
    
    plt.colorbar(ct)
    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    plt.gca().set_aspect('equal')
    
    plt.savefig("scatter_PSIp_JET_95128.pdf", format="pdf", bbox_inches="tight")
    plt.show()  

if plot_evoCon_DiMES==1:
    fig,ax=plt.subplots()
    
    plotcombined=1
    
    if plotcombined==0:
        ax.plot(time[1:], (deconfined[1:]-deconfined[1])/confined[1],'r-',lw=2)
        ax.plot(time[1:], (DiMESdeconfined[1:]-DiMESdeconfined[1])/confined[1],'r:',lw=2)
        
        ax.legend(['total','to DiMES'])
    else:
        ax.plot(time18[1:], (deconfined18[1:]-deconfined18[1])/confined18[1],'r-',lw=2)
        ax.plot(time18a[1:], (deconfined18a[1:]-deconfined18a[1])/confined18a[1],'b-',lw=2)
        ax.plot(time18b[1:], (deconfined18b[1:]-deconfined18b[1])/confined18b[1],'g-',lw=2)
        
        ax.plot(time18[1:], (DiMESdeconfined18[1:]-DiMESdeconfined18[1])/confined18[1],'r:',lw=2)
        ax.plot(time18a[1:], (DiMESdeconfined18a[1:]-DiMESdeconfined18a[1])/confined18a[1],'b:',lw=2)
        ax.plot(time18b[1:], (DiMESdeconfined18b[1:]-DiMESdeconfined18b[1])/confined18b[1],'g:',lw=2)
    
        ax.legend(['A=1','A=10','A=100'])
    
    ax.set(xlabel='$t (\\mathrm{s})$', ylabel='Fraction deconfined')
    ax.grid()
    ax.set_yscale('log')
    
    if plotcombined==0:
        plt.savefig("EvoCon_DIIID_"+run_directory[0]+".png", format="png", bbox_inches="tight")
    else:
        plt.savefig("EvoCon_DIIID_combined.png", format="png", bbox_inches="tight")
    plt.show()   

if plot_histRZ_ext==1:
    
    rbin=np.linspace(min(Rg),max(Rg),num=50)
    zbin=np.linspace(min(Zg),max(Zg),num=50)
    
    rrbin,zzbin=np.meshgrid(rbin,zbin,indexing='ij')
    
    
    plotden=0
    plotconfined=1
    
    if plotden==0:
        if plotconfined==1:
            H, xedges, yedges = np.histogram2d(R[timeind_p,flagCon[timeind_p,:]>0], 
                                          zz[timeind_p,flagCon[timeind_p,:]>0], bins=(rbin, zbin))
        else:
            H, xedges, yedges = np.histogram2d(R[timeind_p,flagCon[timeind_p,:]<1], 
                                          zz[timeind_p,flagCon[timeind_p,:]<1], bins=(rbin, zbin))
    else:
        if plotconfined==1:
            H, xedges, yedges = np.histogram2d(R[timeind_p,flagCon[timeind_p,:]>0], 
                                          zz[timeind_p,flagCon[timeind_p,:]>0], bins=(rbin, zbin),density=True)
        else:
            H, xedges, yedges = np.histogram2d(R[timeind_p,flagCon[timeind_p,:]<1], 
                                          zz[timeind_p,flagCon[timeind_p,:]<1], bins=(rbin, zbin),density=True)
    
    
    fig,ax=plt.subplots()
    
    cmap = plt.cm.viridis

    # Create a modified colormap with 'white' as the "under" color
    cmap_modified = cmap.with_extremes(under='white')
        

    if plotden==1:
        ct=ax.pcolormesh(rrbin, zzbin, H, cmap=cmap_modified,vmin=0.99*np.min(H[H>0]))
    else:
        ct=ax.pcolormesh(rrbin, zzbin, H, cmap=cmap_modified,vmin=1)

    
    #ct=ax.pcolormesh(rrbin, zzbin, H)
    ax.contour(Rg,Zg,-PSIp,25,colors='black',linewidths=0.5)
    ax.contour(Rg,Zg,FLAG,[.5],colors='black')

    
    cb=plt.colorbar(ct)
    
    if plotden==0:
        cb.set_label('$N_{RE}$', rotation=90)
    else:
        cb.set_label('$f_{RE}$', rotation=90)
    
    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    plt.gca().set_aspect('equal')
    
    plt.savefig("RZhist_DIIID.png", format="png", bbox_inches="tight")
    plt.show()
    

if plot_psip==1:
    fig,ax=plt.subplots()
    
    nemin=np.min(-PSIp[:,timeind_g,:])
    nemax=np.max(-PSIp[:,timeind_g,:])
    
    
    ct=ax.contourf(Rg,Zg,-PSIp[:,timeind_g,:],100,cmap=plt.cm.viridis,vmin=nemin,vmax=nemax)
    ax.contour(Rg,Zg,-PSIp[:,timeind_g,:],25,colors='black',linewidths=0.5)
    ax.contour(Rg,Zg,FLAG[:,timeind_g,:],[.5],colors='black')
    
    sc=ax.scatter(R,zz,c=-psiP,s=5,cmap=ct.cmap,vmin=nemin,vmax=nemax,edgecolor='black',linewidth=0.25)
    
    plt.colorbar(ct)
    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    plt.gca().set_aspect('equal')
    
    plt.savefig("scatter_PSIp_JET_95128.pdf", format="pdf", bbox_inches="tight")
    plt.show()  

if plot_ephi==1:
    fig,ax=plt.subplots()
    
    #nemin=np.min(EPHI_g[:,timeind_g,:])
    #nemax=np.max(EPHI_g[:,timeind_g,:])
    
    nemin=np.min(ePHI[timeind_p,flagRE[timeind_p,:]>0])
    nemax=np.max(ePHI[timeind_p,flagRE[timeind_p,:]>0])
    
    #if abs(nemin)>nemax:
    #    nemax=abs(nemin)
    #else:
    #    nemin=-nemax
    
    #ct=ax.contourf(Rg,Zg,EPHI_g[:,timeind_g,:],100,cmap=plt.cm.bwr,vmin=nemin,vmax=nemax)
    #ax.contour(Rg,Zg,PSIp[:,timeind_g,:],25,colors='black',linewidths=0.5)
    #ax.contour(Rg,Zg,FLAG[:,timeind_g,:],[.5],colors='black')
    
    sc=ax.scatter(R[timeind_p,flagRE[timeind_p,:]>0],zz[timeind_p,flagRE[timeind_p,:]>0],
                  c=ePHI[timeind_p,flagRE[timeind_p,:]>0],s=5,edgecolor='black',linewidth=0.25)
    
    plt.colorbar(sc)
    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    plt.gca().set_aspect('equal')
    ax.set(title='$E_\\phi$')
    
    
    plt.savefig("scatter_EPHI_JET_95128.pdf", format="pdf", bbox_inches="tight")
    plt.show()  

if plot_bz==1:
    fig,ax=plt.subplots()
    
    nemin=np.min(bZ[timeind_p,flagRE[timeind_p,:]>0])
    nemax=np.max(bZ[timeind_p,flagRE[timeind_p,:]>0])
    
    if abs(nemin)>nemax:
        nemax=abs(nemin)
    else:
        nemin=-nemax
    
    #ct=ax.contourf(Rg,Zg,EPHI_g[:,timeind_g,:],100,cmap=plt.cm.bwr,vmin=nemin,vmax=nemax)
    #ax.contour(Rg,Zg,PSIp[:,timeind_g,:],25,colors='black',linewidths=0.5)
    #ax.contour(Rg,Zg,FLAG[:,timeind_g,:],[.5],colors='black')
    
    sc=ax.scatter(R[timeind_p,flagRE[timeind_p,:]>0],zz[timeind_p,flagRE[timeind_p,:]>0],
                  c=bZ[timeind_p,flagRE[timeind_p,:]>0],cmap=plt.cm.bwr,
                  s=5,edgecolor='black',linewidth=0.25,vmin=nemin,vmax=nemax)
    
    plt.colorbar(sc)
    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    plt.gca().set_aspect('equal')
    ax.set(title='$B_Z$')
    
    #plt.savefig("scatter_bz_JET_95128.pdf", format="pdf", bbox_inches="tight")
    plt.show()  

if plot_bphi==1:
    fig,ax=plt.subplots()
    
    nemin=np.min(bPHI[timeind_p,flagRE[timeind_p,:]>0])
    nemax=np.max(bPHI[timeind_p,flagRE[timeind_p,:]>0])
    
    if abs(nemin)>nemax:
        nemax=abs(nemin)
    else:
        nemin=-nemax
    
    #ct=ax.contourf(Rg,Zg,EPHI_g[:,timeind_g,:],100,cmap=plt.cm.bwr,vmin=nemin,vmax=nemax)
    #ax.contour(Rg,Zg,PSIp[:,timeind_g,:],25,colors='black',linewidths=0.5)
    #ax.contour(Rg,Zg,FLAG[:,timeind_g,:],[.5],colors='black')
    
    sc=ax.scatter(R[timeind_p,flagRE[timeind_p,:]>0],zz[timeind_p,flagRE[timeind_p,:]>0],
                  c=bPHI[timeind_p,flagRE[timeind_p,:]>0],s=5,edgecolor='black',linewidth=0.25)
    
    plt.colorbar(sc)
    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    plt.gca().set_aspect('equal')
    ax.set(title='$B_\\phi$')
    
    #plt.savefig("scatter_bphi_JET_95128.pdf", format="pdf", bbox_inches="tight")
    plt.show()  

if plot_br==1:
    fig,ax=plt.subplots()
    
    nemin=np.min(bR[timeind_p,(~np.isnan(bR[timeind_p,:]))])
    nemax=np.max(bR[timeind_p,(~np.isnan(bR[timeind_p,:]))])
    
    if abs(nemin)>nemax:
        nemax=abs(nemin)
    else:
        nemin=-nemax
    
    #ct=ax.contourf(Rg,Zg,BR_g[:,timeind_g,:],100,cmap=plt.cm.bwr,vmin=nemin,vmax=nemax)
    #ax.contour(Rg,Zg,BR_g[:,timeind_g,:],25,colors='black',linewidths=0.5)
    #ax.contour(Rg,Zg,PSIp[:,timeind_g,:],25,colors='black',linewidths=0.5)
    #ax.contour(Rg,Zg,FLAG,[.5],colors='black')
    
    
    
    sc=ax.scatter(R[timeind_p,:],zz[timeind_p,:],c=bR[timeind_p,:],s=5,
                  cmap=plt.cm.bwr,vmin=nemin,vmax=nemax,edgecolor='black',linewidth=0.5)
    
    sc=ax.scatter(R[timeind_p,flagRE[timeind_p,:]>0],zz[timeind_p,flagRE[timeind_p,:]>0],
                  c=bR[timeind_p,flagRE[timeind_p,:]>0],cmap=plt.cm.bwr,
                  s=5,edgecolor='black',linewidth=0.25,vmin=nemin,vmax=nemax)
    
    plt.colorbar(sc)
    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    plt.gca().set_aspect('equal')
    ax.set(title='$B_R$')
    
    
    #plt.savefig("scatter_br_JET_95128.pdf", format="pdf", bbox_inches="tight")
    plt.show()  

if plot_ne==1:
    fig,ax=plt.subplots()
    
    nemin=np.min(ne_g[:,timeind_g,:])
    nemax=np.max(ne_g[:,timeind_g,:])
    
    ct=ax.contourf(Rg,Zg,ne_g[:,timeind_g,:],100,cmap=plt.cm.viridis,vmin=nemin,vmax=nemax)
    ax.contour(Rg,Zg,PSIp[:,timeind_g,:],25,colors='black',linewidths=0.5)
    ax.contour(Rg,Zg,FLAG[:,timeind_g,:],[.5],colors='black')
    
    sc=ax.scatter(R[timeind_p,:],zz[timeind_p,:],c=ne[timeind_p,:],s=5,cmap=ct.cmap,vmin=nemin,vmax=nemax,edgecolor='black',linewidth=0.25)
    
    plt.colorbar(ct)
    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    plt.gca().set_aspect('equal')
    
    plt.savefig("scatter_ne_JET_95128.pdf", format="pdf", bbox_inches="tight")
    plt.show()  

if plot_Te==1:
    fig,ax=plt.subplots()
    
    nemin=np.min(Te_g[:,timeind_g,:])
    nemax=np.max(Te_g[:,timeind_g,:])
    
    ct=ax.contourf(Rg,Zg,Te_g[:,timeind_g,:],100,cmap=plt.cm.viridis,vmin=nemin,vmax=nemax)
    ax.contour(Rg,Zg,PSIp[:,timeind_g,:],25,colors='black',linewidths=0.5)
    ax.contour(Rg,Zg,FLAG[:,timeind_g,:],[.5],colors='black')
    
    sc=ax.scatter(R[timeind_p,:],zz[timeind_p,:],c=Te[timeind_p,:],s=5,cmap=ct.cmap,vmin=nemin,vmax=nemax,edgecolor='black',linewidth=0.25)
    
    plt.colorbar(ct)
    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    ax.grid()
    plt.gca().set_aspect('equal')
    
    plt.savefig("scatter_Te_JET_95128.pdf", format="pdf", bbox_inches="tight")
    plt.show()      

if plotyorbit==1:
    fig,ax=plt.subplots(2)
    ax[0].plot(time,yy,color='b')


    ax[0].set(xlabel='t (s)', ylabel='y (m)',
          title='y orbit')
    
    ax[0].legend(['KORC','Analytic'])
    
    ax[0].grid()
    
    ax[1].plot(time,zz,color='b')

    
    ax[1].set(xlabel='t (s)', ylabel='z (m)',
          title='y orbit')
    
    ax[1].grid()
    
    plt.show()