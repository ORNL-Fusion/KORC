import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import h5py
import shapely.geometry as shape
from contourpy import contour_generator
import sys

sys.path.insert(0, '/home/21b/Python/TokDesigner/kernel/divertor')

import EFITutils as EFITutils

qe=1.6022e-19
me=9.1094e-31
c=2.9979e8
mu0=4*np.pi*10**(-7)

#For KORC_RUNS/LOCAL/TEST6 for DIII-D 200236
Bnorm=2.36507606506348
Lnorm=0.000720699602571918
tnorm=2.40399510841569e-12

#%% open EFIT

#filename='/home/21b/ExperimentalInputs/ARC/MARS/42533M20_geqdsk_t=1.631'
filename='/home/21b/KORC/test/d3d_sol_wall/g131041.2403'

geq = EFITutils.geq(filename)
geq.setup_bfield()
#geq.read_gfile(geq.gfile_name)

Zxpt=geq.xpt_info['xpt1']['zx']

R0=geq.g['rmaxis']
Z0=geq.g['zmaxis']
Rlim=geq.g['rlim']
Zlim=geq.g['zlim']

Lim=[Rlim,Zlim]
Lim=np.transpose(Lim)
Lim=shape.polygon.Polygon(Lim)

Rcentr=geq.g['rzero']
Bcentr=geq.g['bcentr']

RRe=geq.g['R']
ZZe=geq.g['Z']
PSIP=geq.g['psirz']

sface=geq.g['qpsi']

ip_sign=geq.g['ip_sign']

PSIP0e=geq.g['ssimag']
PSIPsole=geq.g['ssibry']

wall_file='/home/21b/M3DC1/unstructured/device_data/ITER/first_wall.pts'
with open(wall_file, 'r') as f: lines = f.readlines()
Npts=int(lines[0].strip())

wall=np.zeros([Npts,2])
Rwall=np.zeros([Npts,1])
Zwall=np.zeros([Npts,1])
for ii in range(0,Npts):
    line = lines[ii+1].strip().split()
    Rwall[ii]=float(line[0])
    Zwall[ii]=float(line[1])
    wall[ii,0]=float(line[0])
    wall[ii,1]=float(line[1])
#wall=np.transpose(wall)
wall=shape.polygon.Polygon(wall)


#%% open Yueqiang's file with MARS fields

openMARS=0
if openMARS==1:

    NRm=101
    NZm=161
    NPHIm=6
    
    filename=['']*NPHIm
    
    filename[0]='/home/21b/ExperimentalInputs/ARC/MARS/Root1/BPLASMA_RECTRZ.OUT'
    filename[1]='/home/21b/ExperimentalInputs/ARC/MARS/Root2/BPLASMA_RECTRZ.OUT'
    filename[2]='/home/21b/ExperimentalInputs/ARC/MARS/Root3/BPLASMA_RECTRZ.OUT'
    filename[3]='/home/21b/ExperimentalInputs/ARC/MARS/Root4/BPLASMA_RECTRZ.OUT'
    filename[4]='/home/21b/ExperimentalInputs/ARC/MARS/Root5/BPLASMA_RECTRZ.OUT'
    filename[5]='/home/21b/ExperimentalInputs/ARC/MARS/Root6/BPLASMA_RECTRZ.OUT'
    
    
    RRm=np.zeros([NRm,NZm])
    ZZm=np.zeros([NRm,NZm])
    ReBR=np.zeros([NRm,NPHIm,NZm])
    ImBR=np.zeros([NRm,NPHIm,NZm])
    ReBZ=np.zeros([NRm,NPHIm,NZm])
    ImBZ=np.zeros([NRm,NPHIm,NZm])
    ReBPHI=np.zeros([NRm,NPHIm,NZm])
    ImBPHI=np.zeros([NRm,NPHIm,NZm])
    flagm=np.zeros([NRm,NZm])
    
    for kk in range(0,NPHIm):
        
        with open(filename[kk], 'r') as f: lines = f.readlines()
        
        for ii in range(0,NRm):
            for jj in range(0,NZm):
                
                line = lines[ii*NZm+jj].strip().split()
    
                ReBR[ii,kk,jj]=float(line[2])
                ImBR[ii,kk,jj]=float(line[3])
                ReBZ[ii,kk,jj]=float(line[4])
                ImBZ[ii,kk,jj]=float(line[5])
                ReBPHI[ii,kk,jj]=float(line[6])
                ImBPHI[ii,kk,jj]=float(line[7])
            
                if kk==0:
                    RRm[ii,jj]=float(line[0])
                    ZZm[ii,jj]=float(line[1])
                    point=shape.Point([line[0],line[1]])
                    flagm[ii,jj]=int(Lim.contains(point))
    
    filename[0]='/home/21b/ExperimentalInputs/ARC/MARS/Root1/EPLASMA_RECTRZ.OUT'
    filename[1]='/home/21b/ExperimentalInputs/ARC/MARS/Root2/EPLASMA_RECTRZ.OUT'
    filename[2]='/home/21b/ExperimentalInputs/ARC/MARS/Root3/EPLASMA_RECTRZ.OUT'
    filename[3]='/home/21b/ExperimentalInputs/ARC/MARS/Root4/EPLASMA_RECTRZ.OUT'
    filename[4]='/home/21b/ExperimentalInputs/ARC/MARS/Root5/EPLASMA_RECTRZ.OUT'
    filename[5]='/home/21b/ExperimentalInputs/ARC/MARS/Root6/EPLASMA_RECTRZ.OUT'
    
    ReER=np.zeros([NRm,NPHIm,NZm])
    ImER=np.zeros([NRm,NPHIm,NZm])
    ReEZ=np.zeros([NRm,NPHIm,NZm])
    ImEZ=np.zeros([NRm,NPHIm,NZm])
    ReEPHI=np.zeros([NRm,NPHIm,NZm])
    ImEPHI=np.zeros([NRm,NPHIm,NZm])
    
    for kk in range(0,NPHIm):
        
        with open(filename[kk], 'r') as f: lines = f.readlines()
        
        for ii in range(0,NRm):
            for jj in range(0,NZm):
                
                line = lines[jj*NRm+ii].strip().split()
    
                ReER[ii,kk,jj]=float(line[2])
                ImER[ii,kk,jj]=float(line[3])
                ReEZ[ii,kk,jj]=float(line[4])
                ImEZ[ii,kk,jj]=float(line[5])
                ReEPHI[ii,kk,jj]=float(line[6])
                ImEPHI[ii,kk,jj]=float(line[7])
    
            
    #%% analysis of MARS-F fields
          
    tauA=5.2057e-7
        
    
    #for case-A
    AMP1=4.8572E-05
    phi1=0
    nn1=1
    gamma1=4.38790E-1/tauA
    omega1=4.39625E-1/tauA
    
    AMP2=3.1992E-05
    phi2=0
    nn2=1
    gamma2=2.40979E-1/tauA
    omega2=-5.93756E-1/tauA
    
    AMP3=9.2756E-05
    phi3=0
    nn3=1
    gamma3=1.60252E-1/tauA
    omega3=-1.25379E+0/tauA
    
    AMP4=4.8572E-05
    phi4=0
    nn4=1
    gamma4=1.47962E-1/tauA
    omega4=-1.15931E+0/tauA
    
    AMP5=3.1992E-05
    phi5=0
    nn5=1
    gamma5=4.93633E-2/tauA
    omega5=-1.52892E+0/tauA
    
    AMP6=9.2756E-05
    phi6=0
    nn6=1
    gamma6=1.21977E-2/tauA
    omega6=-2.13030E+0/tauA
        
        
    PHIm=[nn1,nn2,nn3,nn4,nn5,nn6]
    
    AMP=[AMP1,AMP2,AMP3,AMP4,AMP5,AMP6]
    GR=[gamma1,gamma2,gamma3,gamma4,gamma5,gamma6]
    FR=[omega1,omega2,omega3,omega4,omega5,omega6]
    
    
    
    print_points=0
        
    if print_points==1:
        
        B11=np.zeros([NRm,NPHIm,NZm])
        BZ11=np.zeros([NRm,NPHIm,NZm])
    
        for kk in range(0,NPHIm):
            B11[:,kk,:]=AMP[kk]*np.sqrt(ReBR[:,kk,:]**2+ImBR[:,kk,:]**2+ReBPHI[:,kk,:]**2+ImBPHI[:,kk,:]**2+ReBZ[:,kk,:]**2+ImBZ[:,kk,:]**2)
            BZ11[:,kk,:]=AMP[kk]*np.sqrt(ReBZ[:,kk,:]**2+ImBZ[:,kk,:]**2)
            
        tt=0    
            
        Btot=B11[:,0,:]*np.exp(gamma1*tt)+B11[:,1,:]*np.exp(gamma1*tt)+B11[:,2,:]*np.exp(gamma1*tt)
    
        #initial t=0
        #Rq=1.6
        #PHIq=0.
        #Zq=-0.75
    
    
        #cpu
        #Rq=2693.370804461797*Lnorm
        #PHIq=0.4050619006494515
        #Zq=-1477.963615062793*Lnorm
    
        #gpu
        Rq=2644.825260775237*Lnorm
        #PHIq=0.9975833090971448
        PHIq=0
        Zq=-1505.524334132253*Lnorm
        
        
        print('R,PHI,Z')
        print(Rq/Lnorm)
        print(PHIq)
        print(Zq/Lnorm)
    
    
        ReBR1interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ReBR[:,0,:],method='cubic',bounds_error=False, fill_value=None)
        ImBR1interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ImBR[:,0,:],method='cubic',bounds_error=False, fill_value=None)
        ReBPHI1interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ReBPHI[:,0,:],method='cubic',bounds_error=False, fill_value=None)
        ImBPHI1interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ImBPHI[:,0,:],method='cubic',bounds_error=False, fill_value=None)
        ReBZ1interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ReBZ[:,0,:],method='cubic',bounds_error=False, fill_value=None)
        ImBZ1interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ImBZ[:,0,:],method='cubic',bounds_error=False, fill_value=None)
        
        ReBR2interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ReBR[:,1,:],method='cubic',bounds_error=False, fill_value=None)
        ImBR2interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ImBR[:,1,:],method='cubic',bounds_error=False, fill_value=None)
        ReBPHI2interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ReBPHI[:,1,:],method='cubic',bounds_error=False, fill_value=None)
        ImBPHI2interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ImBPHI[:,1,:],method='cubic',bounds_error=False, fill_value=None)
        ReBZ2interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ReBZ[:,1,:],method='cubic',bounds_error=False, fill_value=None)
        ImBZ2interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ImBZ[:,1,:],method='cubic',bounds_error=False, fill_value=None)
        
        ReBR3interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ReBR[:,2,:],method='cubic',bounds_error=False, fill_value=None)
        ImBR3interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ImBR[:,2,:],method='cubic',bounds_error=False, fill_value=None)
        ReBPHI3interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ReBPHI[:,2,:],method='cubic',bounds_error=False, fill_value=None)
        ImBPHI3interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ImBPHI[:,2,:],method='cubic',bounds_error=False, fill_value=None)
        ReBZ3interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ReBZ[:,2,:],method='cubic',bounds_error=False, fill_value=None)
        ImBZ3interp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),ImBZ[:,2,:],method='cubic',bounds_error=False, fill_value=None)
        
        ReBR1q=ReBR1interp((Rq,Zq))
        ImBR1q=ImBR1interp((Rq,Zq))
        ReBPHI1q=ReBPHI1interp((Rq,Zq))
        ImBPHI1q=ImBPHI1interp((Rq,Zq))
        ReBZ1q=ReBZ1interp((Rq,Zq))
        ImBZ1q=ImBZ1interp((Rq,Zq))
    
        ReBR2q=ReBR2interp((Rq,Zq))
        ImBR2q=ImBR2interp((Rq,Zq))
        ReBPHI2q=ReBPHI2interp((Rq,Zq))
        ImBPHI2q=ImBPHI2interp((Rq,Zq))
        ReBZ2q=ReBZ2interp((Rq,Zq))
        ImBZ2q=ImBZ2interp((Rq,Zq))
        
        ReBR3q=ReBR3interp((Rq,Zq))
        ImBR3q=ImBR3interp((Rq,Zq))
        ReBPHI3q=ReBPHI3interp((Rq,Zq))
        ImBPHI3q=ImBPHI3interp((Rq,Zq))
        ReBZ3q=ReBZ3interp((Rq,Zq))
        ImBZ3q=ImBZ3interp((Rq,Zq))
    
        print('complex B1')
        print(ReBR1q/Bnorm)
        print(ImBR1q/Bnorm)
        print(ReBPHI1q/Bnorm)
        print(ImBPHI1q/Bnorm)
        print(ReBZ1q/Bnorm)
        print(ImBZ1q/Bnorm)
        
        print('complex B2')
        print(ReBR2q/Bnorm)
        print(ImBR2q/Bnorm)
        print(ReBPHI2q/Bnorm)
        print(ImBPHI2q/Bnorm)
        print(ReBZ2q/Bnorm)
        print(ImBZ2q/Bnorm)
    
        cshift1=np.cos(PHIq*PHIm[0])
        sshift1=np.sin(PHIq*PHIm[0])
        cshift2=np.cos(PHIq*PHIm[1])
        sshift2=np.sin(PHIq*PHIm[1])
        cshift3=np.cos(PHIq*PHIm[2])
        sshift3=np.sin(PHIq*PHIm[2])
    
        print('toroidal shift')
        print(cshift1)
        print(sshift1)
    
        amp_scale=1
        bsat=100
        #time=2.0999996893030134E-005 #cpu
        #time=2.0999993786060270E-005 #gpu
        time=4e-5
    
        print('time evo')
        print(amp_scale*AMP[0]*bsat/(1+bsat*np.exp(-GR[0]*time)))
    
        B11R=amp_scale*AMP[0]*(ReBR1q*cshift1-ImBR1q*sshift1)*bsat/(1+bsat*np.exp(-GR[0]*time))
        B12R=amp_scale*AMP[1]*(ReBR2q*cshift2-ImBR2q*sshift2)*bsat/(1+bsat*np.exp(-GR[1]*time))
        B13R=amp_scale*AMP[2]*(ReBR3q*cshift3-ImBR3q*sshift3)*bsat/(1+bsat*np.exp(-GR[2]*time))
        
        B11PHI=amp_scale*AMP[0]*(ReBPHI1q*cshift1-ImBPHI1q*sshift1)*bsat/(1+bsat*np.exp(-GR[0]*time))
        B12PHI=amp_scale*AMP[1]*(ReBPHI2q*cshift2-ImBPHI2q*sshift2)*bsat/(1+bsat*np.exp(-GR[1]*time))
        B13PHI=amp_scale*AMP[2]*(ReBPHI3q*cshift3-ImBPHI3q*sshift3)*bsat/(1+bsat*np.exp(-GR[2]*time))
        
        B11Z=amp_scale*AMP[0]*(ReBZ1q*cshift1-ImBZ1q*sshift1)*bsat/(1+bsat*np.exp(-GR[0]*time))
        B12Z=amp_scale*AMP[1]*(ReBZ2q*cshift2-ImBZ2q*sshift2)*bsat/(1+bsat*np.exp(-GR[1]*time))
        B13Z=amp_scale*AMP[2]*(ReBZ3q*cshift3-ImBZ3q*sshift3)*bsat/(1+bsat*np.exp(-GR[2]*time))
    
        print('real B1')
        print(B11R/Bnorm)
        print(B11PHI/Bnorm)
        print(B11Z/Bnorm)
        
        print('real B2')
        print(B12R/Bnorm)
        print(B12PHI/Bnorm)
        print(B12Z/Bnorm)
        
        print('real B3')
        print(B13R/Bnorm)
        print(B13PHI/Bnorm)
        print(B13Z/Bnorm)
    
        B11R_2D=amp_scale*AMP[0]*(ReBR[:,0,:]*cshift1-ImBR[:,0,:]*sshift1)*bsat/(1+bsat*np.exp(-GR[0]*time))
        B12R_2D=amp_scale*AMP[1]*(ReBR[:,1,:]*cshift2-ImBR[:,1,:]*sshift2)*bsat/(1+bsat*np.exp(-GR[1]*time))
        B13R_2D=amp_scale*AMP[2]*(ReBR[:,2,:]*cshift3-ImBR[:,2,:]*sshift3)*bsat/(1+bsat*np.exp(-GR[2]*time))
        
        B11PHI_2D=amp_scale*AMP[0]*(ReBPHI[:,0,:]*cshift1-ImBPHI[:,0,:]*sshift1)*bsat/(1+bsat*np.exp(-GR[0]*time))
        B12PHI_2D=amp_scale*AMP[1]*(ReBPHI[:,1,:]*cshift2-ImBPHI[:,1,:]*sshift2)*bsat/(1+bsat*np.exp(-GR[1]*time))
        B13PHI_2D=amp_scale*AMP[2]*(ReBPHI[:,2,:]*cshift3-ImBPHI[:,2,:]*sshift3)*bsat/(1+bsat*np.exp(-GR[2]*time))
        
        B11Z_2D=amp_scale*AMP[0]*(ReBZ[:,0,:]*cshift1-ImBZ[:,0,:]*sshift1)*bsat/(1+bsat*np.exp(-GR[0]*time))
        B12Z_2D=amp_scale*AMP[1]*(ReBZ[:,1,:]*cshift2-ImBZ[:,1,:]*sshift2)*bsat/(1+bsat*np.exp(-GR[1]*time))
        B13Z_2D=amp_scale*AMP[2]*(ReBZ[:,2,:]*cshift3-ImBZ[:,2,:]*sshift3)*bsat/(1+bsat*np.exp(-GR[2]*time))
    

#%% interpolate Eric's fields onto Yueqiang's mesh

PSIPinterp=sci.interpolate.RegularGridInterpolator((RRe,ZZe),PSIP,method='cubic',bounds_error=False, fill_value=None)

if openMARS==1:    
    PSIPm=PSIPinterp((RRm,ZZm))
else:
    PSIPm=PSIP
    ZZm,RRm=np.meshgrid(ZZe,RRe)
    NRm=np.shape(RRm[:,0])[0]
    NZm=np.shape(ZZm[0,:])[0]
    
    flagm=np.zeros([NRm,NZm])
    for ii in range(0,NRm):
        for jj in range(0,NZm):
            point=shape.Point([RRm[ii,jj],ZZm[ii,jj]])
            #flagm[ii,jj]=int(Lim.contains(point))
            flagm[ii,jj]=int(wall.contains(point))
    
PSIP0=np.min(PSIPm)
PSIPlim=PSIPinterp((Rlim,Zlim))
PSIPwall=np.max(PSIPlim)
PSIPlim=np.min(PSIPlim)

#%% analysis of axisymmetric fields

dPSIdR,dPSIdZ=np.gradient(PSIPm,RRm[:,0],ZZm[0,:]);

#psip_conv=-1/(2*np.pi) #for DIII-D
psip_conv=-1 #for ARC

BRm=psip_conv*dPSIdZ/RRm
BZm=-psip_conv*dPSIdR/RRm
BPHIm=np.zeros(np.shape(BRm))
for ii in range(0,NRm):
    BPHIm[ii,:]=Rcentr*Bcentr/RRm[ii]

PSIres=256

PSIParray=np.linspace(PSIP0,PSIPlim,PSIres)
PSINarray=(PSIParray-PSIP0)/(PSIPlim-PSIP0)

PSIT=np.zeros(PSIres)
for jj in range (0,PSIres):
    BPHItmp=BPHIm.copy()
    BPHItmp[flagm<1]=0
    BPHItmp[PSIPm>PSIParray[jj]]=0
    PSIT[jj]=np.trapezoid(RRm[:,0],np.trapezoid(ZZm[0,:],BPHItmp,2))

sfac=np.gradient(PSIT,2*np.pi*PSIParray/psip_conv);

q_debug=0
if q_debug==1:
    fig,ax=plt.subplots()
    ax.plot(PSIParray)
    ax.plot(PSIT)
    
    fig,ax=plt.subplots()
    BPHItmp=BPHIm.copy()
    BPHItmp[flagm<1]=0
    BPHItmp[PSIPm>PSIParray[100]]=0
    ax.contourf(RRm,ZZm,BPHItmp)
    

if print_points==1:
    
    BRminterp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),BRm,method='cubic',bounds_error=False, fill_value=None)
    BPHIminterp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),BPHIm,method='cubic',bounds_error=False, fill_value=None)
    BZminterp=sci.interpolate.RegularGridInterpolator((RRm[:,0],ZZm[0,:]),BZm,method='cubic',bounds_error=False, fill_value=None)
    
    BRmq=BRminterp((Rq,Zq))
    BPHImq=BPHIminterp((Rq,Zq))
    BZmq=BZminterp((Rq,Zq))

    print('B0')
    print(BRmq/Bnorm)
    print(BPHImq/Bnorm)
    print(BZmq/Bnorm)

    BRtot=BRmq+B11R+B12R+B13R
    BPHItot=BPHImq+B11PHI+B12PHI+B13PHI
    BZtot=BZmq+B11Z+B12Z+B13Z
        
    print('total B')
    print(BRtot/Bnorm)
    print(BPHItot/Bnorm)
    print(BZtot/Bnorm)

#%% writing fields to hdf5 file

create_input=1

if create_input==1:
    
    #filename='ARC_42533M20_MARS.h5'
    filename='ITER_IMAS.h5'

    with h5py.File(filename, "w") as f:
        dset=f.create_dataset('NR',(1,),dtype='i')
        dset[0]=NRm
        dset=f.create_dataset('NZ',(1,),dtype='i')
        dset[0]=NZm
        dset=f.create_dataset('R',(NRm,),dtype='f')
        dset[:]=RRm[:,0]
        dset=f.create_dataset('Z',(NZm,),dtype='f')
        dset[:]=ZZm[0,:]
        dset=f.create_dataset('FLAG',(NZm,NRm),dtype='i')
        dset[:,:]=np.transpose(flagm[:,:])
        dset=f.create_dataset('PSIp',(NZm,NRm),dtype='f')
        dset[:,:]=np.transpose(PSIPm[:,:])
        if openMARS==1: 
            dset=f.create_dataset('NPHI',(1,),dtype='i')
            dset[0]=NPHIm
            dset=f.create_dataset('PHI',(NPHIm,),dtype='f')
            dset[:]=PHIm[:]
            dset=f.create_dataset('ReBR',(NZm,NPHIm,NRm),dtype='f')
            dset[:,:]=ReBR.transpose((2,1,0))
            dset=f.create_dataset('ImBR',(NZm,NPHIm,NRm),dtype='f')
            dset[:,:]=ImBR.transpose((2,1,0))
            dset=f.create_dataset('ReBZ',(NZm,NPHIm,NRm),dtype='f')
            dset[:,:]=ReBZ.transpose((2,1,0))
            dset=f.create_dataset('ImBZ',(NZm,NPHIm,NRm),dtype='f')
            dset[:,:]=ImBZ.transpose((2,1,0))
            dset=f.create_dataset('ReBPHI',(NZm,NPHIm,NRm),dtype='f')
            dset[:,:]=ReBPHI.transpose((2,1,0))
            dset=f.create_dataset('ImBPHI',(NZm,NPHIm,NRm),dtype='f')
            dset[:,:]=ImBPHI.transpose((2,1,0))
            dset=f.create_dataset('ReER',(NZm,NPHIm,NRm),dtype='f')
            dset[:,:]=ReER.transpose((2,1,0))
            dset=f.create_dataset('ImER',(NZm,NPHIm,NRm),dtype='f')
            dset[:,:]=ImER.transpose((2,1,0))
            dset=f.create_dataset('ReEZ',(NZm,NPHIm,NRm),dtype='f')
            dset[:,:]=ReEZ.transpose((2,1,0))
            dset=f.create_dataset('ImEZ',(NZm,NPHIm,NRm),dtype='f')
            dset[:,:]=ImEZ.transpose((2,1,0))
            dset=f.create_dataset('ReEPHI',(NZm,NPHIm,NRm),dtype='f')
            dset[:,:]=ReEPHI.transpose((2,1,0))
            dset=f.create_dataset('ImEPHI',(NZm,NPHIm,NRm),dtype='f')
            dset[:,:]=ImEPHI.transpose((2,1,0))
        dset=f.create_dataset('Bo',(1,),dtype='f')
        dset[0]=Bcentr
        dset=f.create_dataset('Ro',(1,),dtype='f')
        dset[0]=Rcentr
        dset=f.create_dataset('Zo',(1,),dtype='f')
        dset[0]=0 
        dset=f.create_dataset('Eo',(1,),dtype='f')
        dset[0]=0
        dset=f.create_dataset('PSIP0',(1,),dtype='f')
        dset[0]=PSIP0
        dset=f.create_dataset('PSIPlim',(1,),dtype='f')
        dset[0]=PSIPlim
        if openMARS==1:  
            dset=f.create_dataset('AMP',(NPHIm,),dtype='f')
            dset[:]=AMP
            dset=f.create_dataset('GR',(NPHIm,),dtype='f')
            dset[:]=GR
            dset=f.create_dataset('FR',(NPHIm,),dtype='f')
            dset[:]=FR
        dset=f.create_dataset('Zxpt',(1,),dtype='f')
        dset[0]=Zxpt
    
#%% Plotting

plt.rcParams['contour.negative_linestyle'] = 'solid'
plt.rcParams['font.size'] = 16  # Default font size for all text
plt.rcParams['axes.titlesize'] = 16  # Font size for axes titles
plt.rcParams['axes.labelsize'] = 16  # Font size for x and y labels
plt.rcParams['xtick.labelsize'] = 16  # Font size for x-axis tick labels
plt.rcParams['ytick.labelsize'] = 16  # Font size for y-axis tick labels
plt.rcParams['legend.fontsize'] = 16  # Font size for legend text
plt.rcParams['figure.titlesize'] = 16 # Font size for figure titles

RZcont=1
RZlogcont=0
tmpcut=0
plotfRE=0
plotGR=0
plotNREscaling=0


if plotNREscaling==1:
    fig,ax=plt.subplots()
    
    NRE=[1,10,100,1000,10000,100000]
    runtime=[7.99,7.97,17.6,17.6,63.3,522]
    
    ax.plot(NRE,runtime,'o-')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.grid()
    ax.set(xlabel='$\\mathrm{N_{RE}}$')
    ax.set(ylabel='$\\mathrm{ runtime\,(s) }$')
    
    #plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

if plotGR==1:
    fig,ax=plt.subplots()

    time=10**(np.linspace(-7,-3,100))
    
    b0=1e-4
    bsat=100
    quas_fac=10**8
    
    btmp=b0*np.exp(GR[0]*time)
    btmp[btmp>bsat*b0]=bsat*b0
    
    btmp1=bsat/(1+bsat*np.exp(-GR[0]*time))*b0
    btmp2=bsat/(1+bsat*np.exp(-GR[1]*time))*b0
    btmp3=bsat/(1+bsat*np.exp(-GR[2]*time))*b0
    
    btmp11=GR[0]/quas_fac/(1+(bsat-1)*np.exp(-GR[0]*time))
    btmp12=GR[1]/quas_fac/(1+(bsat-1)*np.exp(-GR[1]*time))
    btmp13=GR[2]/quas_fac/(1+(bsat-1)*np.exp(-GR[2]*time))
    btmp14=GR[3]/quas_fac/(1+(bsat-1)*np.exp(-GR[3]*time))
    btmp15=GR[4]/quas_fac/(1+(bsat-1)*np.exp(-GR[4]*time))
    btmp16=GR[5]/quas_fac/(1+(bsat-1)*np.exp(-GR[5]*time))
    

    #ax.plot(time,btmp)
    ax.plot(time,btmp11,lw=2)
    ax.plot(time,btmp12,lw=2)
    ax.plot(time,btmp13,lw=2)
    ax.plot(time,btmp14,lw=2)
    ax.plot(time,btmp15,lw=2)
    ax.plot(time,btmp16,lw=2)
    
    ax.axis([0,4e-4,1e-6,1e-2])
    
    #ax.set_xscale('log')
    ax.set_yscale('log')
    
    #ax.legend(['n=1, Root 1','n=1, Root 2','n=2, Root 1'])
    
    ax.grid()
    ax.set(xlabel='$t (\\mathrm{s})$')
    
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
    plt.savefig("B1_amp_evo.png", format="png", bbox_inches="tight")
    plt.show()


if tmpcut==1:
    fig,ax=plt.subplots()
    
    #tmpfld=RRm
    #tmpfld=BPHIm
    tmpfld=sfac
    
    #ax.plot(RRm)
    #ax.plot(RRe,tmpfld[33,:])
    ax.plot(PSINarray,sfac)
    ax.plot(PSINarray,sface)
    
    #R0tmp=RRe[20]
    #B0tmp=BPHIe[33,20]
    
    #ax.plot(RRe,R0tmp*B0tmp/RRe)
    
    ax.grid()
    #ax.set(xlabel='$R (\\mathrm{m})$')
    
    ax.axis([0,1,-0,10])
    
    plt.show()

if RZlogcont==1:
    fig,ax=plt.subplots()
    
    #tmpfld=np.sqrt(ReBR*ReBR+ImBR*ImBR)
    tmpfld=BRm
    #tmpfld=PSIPm
    #tmpfld=B11[:,2,:]
    #tmpfld=Btot
    
    
    tmpfld[flagm<1]=0
    
    levs10=np.linspace(np.log10(np.min(tmpfld[flagm>0])),np.log10(np.max(tmpfld[flagm>0])),50)
    levs=10**levs10
    
    #ct=ax.contourf(RReg,ZZeg,tmpfld,vmin=np.min(tmpfld),vmax=np.max(tmpfld))
    ct=ax.contourf(RRm,ZZm,np.log10(tmpfld),levels=levs10,
                   vmin=np.min(np.log10(tmpfld[flagm>0])),vmax=np.max(np.log10(tmpfld[flagm>0])))
    ax.plot(Rlim,Zlim,'w-')
    
    plt.colorbar(ct)
    ax.set(xlabel='R', ylabel='Z')
    
    ax.grid()
    plt.gca().set_aspect('equal')
    
    plt.show()

if RZcont==1:
    fig,ax=plt.subplots()
    
    #lo=0.0007206996025719176
    #bo=2.3650760650634766
    
    #tmpfld=np.sqrt(ReBR*ReBR+ImBR*ImBR)
    #tmpfld=BPHIm.copy()
    #tmpfld=B13R_2D.copy()
    #tmpfld=B11[:,0,:]
    #tmpfld=dPSIdR
    #tmpfld=dPSIdZ
    tmpfld=PSIPm.copy()
    #tmpfld=dPSIdZ.copy()
    #tmpfld=nREm
    #tmpfld=ReER.copy()
    #tmpfld=BZm.copy()

    #tmpfld=tmpfld[:,0,:]

    #tmpfld=tmpfld/(bo*lo**2)
    
    tmpfld=np.ma.masked_where(flagm<1,tmpfld)
    
    tmpfld1=PSIPm.copy()
    
    vmin=np.min(tmpfld)
    vmax=np.max(tmpfld)
    #vmax=-1183954 #debug PSIP
    #vmin=-261 #debug dPSIPdR
    #vmax=261
    #vmin=-4384 #debug dPSIPdZ
    #vmax=4384
    
    if vmin*vmax<0:
        if vmax>np.abs(vmin):
            vmin=-vmax
        else:
            vmax=np.abs(vmin)
    
    
    levs=np.linspace(vmin,vmax,50)
    #levs=10**levs10
    
    #ct=ax.contourf(RReg,ZZeg,tmpfld,vmin=np.min(tmpfld),vmax=np.max(tmpfld))
    if vmin*vmax<0:
        ct=ax.contourf(RRm,ZZm,tmpfld,levs,cmap=plt.cm.bwr)
    else:
        ct=ax.contourf(RRm,ZZm,tmpfld,levs)
        
    ax.contour(RRm,ZZm,tmpfld1,[PSIPsole])
        
    #ax.contour(RRm,ZZm,tmpfld1,20,colors='black',linewidths=1)
    
    #contour=contour_generator(RRm,ZZm,RHON)
    #RZcontour=contour.lines(1)[0]
    
    #ax.plot([1.6],[-0.75],'rx')
    #ax.plot(R0,Z0,'k+')
    
    #ax.plot(RZcontour[:,0],RZcontour[:,1],'k--')
    #ax.plot(Rlim,Zlim,'k-')
    ax.plot(Rwall,Zwall,'k-')
    
    cb=plt.colorbar(ct)
    cb.formatter.set_powerlimits((0, 0))
    cb.formatter.set_useMathText(True)


    ax.set(xlabel='$R (\\mathrm{m})$', ylabel='$Z (\\mathrm{m})$')
    
    ax.grid()
    plt.gca().set_aspect('equal')
    
    plt.savefig("B0Z_evo.png", format="png", bbox_inches="tight")
    plt.show()