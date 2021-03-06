% getderivedNIMRODFields.m

%% Initialization

clear all
format longE

% load desired matlab datafile generated by getNIMRODFields.m
runnum=00000;
filename=strcat('NIMROD_MST',num2str(runnum,'%05i'),'.mat');
load(filename)

%% Manipulation

% making matricies periodic in phi

S.PHI=[S.PHI 2*pi];

BRtmp=zeros(size(S.R,2),size(S.PHI,2),size(S.Z,2));
BPHItmp=zeros(size(S.R,2),size(S.PHI,2),size(S.Z,2));
BZtmp=zeros(size(S.R,2),size(S.PHI,2),size(S.Z,2));
VRtmp=zeros(size(S.R,2),size(S.PHI,2),size(S.Z,2));
VPHItmp=zeros(size(S.R,2),size(S.PHI,2),size(S.Z,2));
VZtmp=zeros(size(S.R,2),size(S.PHI,2),size(S.Z,2));
JRtmp=zeros(size(S.R,2),size(S.PHI,2),size(S.Z,2));
JPHItmp=zeros(size(S.R,2),size(S.PHI,2),size(S.Z,2));
JZtmp=zeros(size(S.R,2),size(S.PHI,2),size(S.Z,2));
DIFFtmp=zeros(size(S.R,2),size(S.PHI,2),size(S.Z,2));
FLAGtmp=zeros(size(S.R,2),size(S.PHI,2),size(S.Z,2));

BRtmp(:,1:end-1,:)=S.BR;
BPHItmp(:,1:end-1,:)=S.BPHI;
BZtmp(:,1:end-1,:)=S.BZ;
VRtmp(:,1:end-1,:)=S.VR;
VPHItmp(:,1:end-1,:)=S.VPHI;
VZtmp(:,1:end-1,:)=S.VZ;
JRtmp(:,1:end-1,:)=S.JR;
JPHItmp(:,1:end-1,:)=S.JPHI;
JZtmp(:,1:end-1,:)=S.JZ;
DIFFtmp(:,1:end-1,:)=S.DIFFPROF;
FLAGtmp(:,1:end-1,:)=S.FLAG;

BRtmp(:,end,:)=BRtmp(:,1,:);
BPHItmp(:,end,:)=BPHItmp(:,1,:);
BZtmp(:,end,:)=BZtmp(:,1,:);
VRtmp(:,end,:)=VRtmp(:,1,:);
VPHItmp(:,end,:)=VPHItmp(:,1,:);
VZtmp(:,end,:)=VZtmp(:,1,:);
JRtmp(:,end,:)=JRtmp(:,1,:);
JPHItmp(:,end,:)=JPHItmp(:,1,:);
JZtmp(:,end,:)=JZtmp(:,1,:);
DIFFtmp(:,end,:)=DIFFtmp(:,1,:);
FLAGtmp(:,end,:)=FLAGtmp(:,1,:);

S.BR=BRtmp;
S.BPHI=BPHItmp;
S.BZ=BZtmp;
S.VR=VRtmp;
S.VPHI=VPHItmp;
S.VZ=VZtmp;
S.JR=JRtmp;
S.JPHI=JPHItmp;
S.JZ=JZtmp;
S.DIFFPROF=DIFFtmp;
S.FLAG=FLAGtmp;

mu0=4*pi*10^(-7);

% from nimrod.in input file
elecd=5.453;

eta=elecd*mu0;

S.ER=zeros(size(S.BR));
S.EPHI=zeros(size(S.BR));
S.EZ=zeros(size(S.BR));

S.ER=-(S.VPHI.*S.BZ-S.VZ.*S.BPHI)+eta*S.DIFFPROF.*S.JR;
S.EPHI=-(S.VZ.*S.BR-S.VR.*S.BZ)+eta*S.DIFFPROF.*S.JPHI;
S.EZ=-(S.VR.*S.BPHI-S.VPHI.*S.BR)+eta*S.DIFFPROF.*S.JZ;

JPHI_AS=squeeze(S.JPHI(:,1,:));

IP=trapz(S.Z,trapz(S.R,JPHI_AS));

%% Plotting

plotE=1;

if plotE==1
    fig=figure('position',[100,100,900,250]);

    subplot(1,3,1)
    contourf(S.R,S.Z,squeeze(S.ER(:,1,:))')
    daspect([1,1,1])
    colormap redblue
    colorbar
    maxER=max(S.ER(:));
    minER=min(S.ER(:));
    maxER=max(maxER,abs(minER));
    caxis([-maxER,maxER])

    subplot(1,3,2)
    contourf(S.R,S.Z,squeeze(S.EPHI(:,1,:))')
    daspect([1,1,1])
    colormap redblue
    colorbar;
    maxEPHI=max(S.EPHI(:));
    minEPHI=min(S.EPHI(:));
    maxEPHI=max(maxEPHI,abs(minEPHI));
    caxis([-maxEPHI,maxEPHI])

    subplot(1,3,3)
    contourf(S.R,S.Z,squeeze(S.EZ(:,1,:))')
    daspect([1,1,1])
    colormap redblue
    colorbar;
    maxEZ=max(S.EZ(:));
    minEZ=min(S.EZ(:));
    maxEZ=max(maxEZ,abs(minEZ));
    caxis([-maxEZ,maxEZ])
end

%% Save

ifsave=1;

if ifsave==1
    save(strcat('NIMRODderived_MST',num2str(id,'%05i')))
end

%% Write PSPLINE input file

ifPSPLINE=0;
ifNoE=1;

if ifPSPLINE==1
    if runnum==00000
   
        AS_BR=squeeze(S.BR(:,1,:));
        AS_BPHI=squeeze(S.BPHI(:,1,:));
        AS_BZ=squeeze(S.BZ(:,1,:));
        AS_ER=squeeze(S.ER(:,1,:));
        AS_EPHI=squeeze(S.EPHI(:,1,:));
        AS_EZ=squeeze(S.EZ(:,1,:));
        AS_FLAG=squeeze(S.FLAG(:,1,:));
        
        if ifNoE==0
            
            filename=strcat('NIMROD_MST_',num2str(runnum,'%05i'),'_AS.h5');
            
            delete(filename);
            
            fields2hdf(S.R,[],S.Z,AS_BR,AS_BPHI,AS_BZ,AS_ER,AS_EPHI,AS_EZ,...
                [],AS_FLAG,[],[],[],filename,S.Bo,0,S.Ro,S.Zo); 
        else
            filename=strcat('NIMROD_MST_',num2str(runnum,'%05i'),'_AS_NoE.h5');
            
            delete(filename);
            
            fields2hdf(S.R,[],S.Z,AS_BR,AS_BPHI,AS_BZ,[],[],[],...
                [],AS_FLAG,[],[],[],filename,S.Bo,0,S.Ro,S.Zo);             
        end
    else
        if ifNoE==0
            filename=strcat('NIMROD_MST_',num2str(runnum,'%05i'),'.h5');
            
            delete(filename);
            
            fields2hdf(S.R,S.PHI,S.Z,S.BR,S.BPHI,S.BZ,S.ER,S.EPHI,S.EZ,...
                [],S.FLAG,[],[],[],filename,S.Bo,0,S.Ro,S.Zo); 
        else
            filename=strcat('NIMROD_MST_',num2str(runnum,'%05i'),'_noE.h5');            
            
            delete(filename);
            
            fields2hdf(S.R,S.PHI,S.Z,S.BR,S.BPHI,S.BZ,[],[],[],...
                [],S.FLAG,[],[],[],filename,S.Bo,0,S.Ro,S.Zo);             
        end
    end
end







