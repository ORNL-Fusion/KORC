format longE

%% load data from file

join=0;
if(join==1)
    
    load('../TEST29d/29d/ST.mat')
    ST29d=ST;
    load('../TEST29d/29d1/ST.mat')
    ST29d1=ST;
    
    ST.time=cat(2,ST29d.time,ST29d1.time);

    ST.data.sp1.Y=cat(3,ST29d.data.sp1.Y,ST29d1.data.sp1.Y);
    ST.data.sp1.V=cat(3,ST29d.data.sp1.V,ST29d1.data.sp1.V);
    ST.data.sp1.B=cat(3,ST29d.data.sp1.B,ST29d1.data.sp1.B);
    ST.data.sp1.E=cat(3,ST29d.data.sp1.E,ST29d1.data.sp1.E);
    ST.data.sp1.g=cat(2,ST29d.data.sp1.g,ST29d1.data.sp1.g);
    ST.data.sp1.eta=cat(2,ST29d.data.sp1.eta,ST29d1.data.sp1.eta);
    ST.data.sp1.flag=cat(2,ST29d.data.sp1.flag,ST29d1.data.sp1.flag);
    ST.data.sp1.PSIp=cat(2,ST29d.data.sp1.PSIp,ST29d1.data.sp1.PSIp);
    ST.data.sp1.ne=cat(2,ST29d.data.sp1.ne,ST29d1.data.sp1.ne);

end

%% extract data from ST structure

c=ST.params.scales.v;
me=ST.params.scales.m;
qe=ST.params.scales.q;
e0=8.8542e-12;

time=ST.time;

R=squeeze(ST.data.sp1.Y(:,1,:));
PHI=squeeze(ST.data.sp1.Y(:,2,:));
Z=squeeze(ST.data.sp1.Y(:,3,:));

%R0=ST.params.fields.Ro;
%Z0=ST.params.fields.Zo;
R0=ST.params.species.Ro;
Z0=ST.params.species.Zo;
rm=sqrt((R-R0).^2+(Z-Z0).^2);
theta=atan2((Z-Z0),(R-R0));

BR=squeeze(ST.data.sp1.B(:,1,:));
BPHI=squeeze(ST.data.sp1.B(:,2,:));
BZ=squeeze(ST.data.sp1.B(:,3,:));

EPHI=squeeze(ST.data.sp1.E(:,2,:));

%PSIp=ST.data.sp1.PSIp;

%lam=ST.params.fields.lambda;
%q0=ST.params.fields.qo;
%B0=ST.params.fields.Bo;
%PSIanaly=R*lam^2*B0./(2*q0*(R0+rm.*cos(theta))).*log(1+(rm/lam).^2);

BMAG=sqrt(BR.^2+BPHI.^2+BZ.^2);

ppll=squeeze(ST.data.sp1.V(:,1,:));
mu=squeeze(ST.data.sp1.V(:,2,:));

%momgcv=ppll.*R.*BPHI./BMAG;
%momgcb=qe*PSIp;
%momgc=momgcv+momgcb;

pmag=sqrt(ppll.^2+mu.*BMAG*2*me);

gam=ST.data.sp1.g(:,:);
eta=ST.data.sp1.eta(:,:);
xi=cos(deg2rad(eta));

vpll=ppll./(gam*me);

flagCon=ST.data.sp1.flagCon;
flagCol=ST.data.sp1.flagCol;
%flagRE=ST.data.sp1.flagRE;

Ipart=qe*vpll.*BPHI./BMAG;
Ipart(flagCon==0)=0;
Ipart(flagCol==0)=0;
I=sum(Ipart,1);

gam(flagCon==0)=0;
gam(flagCol==0)=0;
E=me*c^2/qe*sum(gam);
KE=me*c^2/qe*sum((gam-1).*flagCon.*flagCol);

AveE=me*c^2/qe*sum(gam.*flagCon.*flagCol)/size(R,1);
AveE(isnan(AveE))=0;
AveKE=me*c^2/qe*sum((gam-1).*flagCon.*flagCol)/size(R,1);
AveKE(isnan(AveKE))=0;

dgamdt=-qe^4*BMAG.^2/(6*pi*e0*(me*c)^3).*(gam.^2-1).*(1-xi.^2);
dgam=dgamdt*ST.params.simulation.dt;

gamanaly=zeros(size(time));
gamanaly(1)=gam(1,1);
for i=2:size(time,2)
    gamanaly(i)=gamanaly(i-1)+dgam(1,i-1);
end

dpplldt=-ppll.*(1-xi.^2).*(gam-1./gam)./(6*pi*e0*(me*c)^3./(qe^4*BMAG.^2));
dmudt=-2*mu./(6*pi*e0*(me*c)^3./(qe^4*BMAG.^2)).*(gam.*(1-xi.^2)+xi.^2./gam);

%ne=ST.data.sp1.ne;

%% construct background toroidal magnetic vector potential

if (strcmp(ST.params.simulation.field_model,'ANALYTICAL'))
    lam=ST.params.fields.lambda;
    R0=ST.params.fields.Ro;
    Z0=ST.params.fields.Zo; 
    a=ST.params.fields.a;
    B0=ST.params.fields.Bo;
    q0=ST.params.fields.qo;
    
    RF=linspace(R0-a,R0+a,100);
    ZF=linspace(Z0-a,Z0+a,100); 

    rm=zeros(100);
    PSIP=zeros(100);  
    FLAG=ones(100);    
    for ii=1:100
       for jj=1:100
           rm(ii,jj)=sqrt((RF(ii)-R0)^2+(ZF(jj)-Z0)^2);
           PSIP(ii,jj)=lam^2*B0./(2*q0)*log(1+(rm(ii,jj)/lam)^2);
           if rm(ii,jj)>a
               FLAG(ii,jj)=0;
           end
       end
    end
    
elseif (strcmp(ST.params.simulation.field_model{1},'EXTERNAL') &&...
        ST.params.fields.Axisymmetric==0)
    PSIP=ST.params.fields.psi_p;
    FLAG=ST.params.fields.Flag;

    RF=ST.params.fields.R;
    ZF=ST.params.fields.Z;

    RF1=RF(1);
    ZF1=ZF(1);

    skipres=2^0;

    for i=skipres:skipres:size(RF)
        RF1=cat(1,RF1,RF(i));
    end
    for i=skipres:skipres:size(ZF)
        ZF1=cat(1,ZF1,ZF(i));
    end
    
elseif (strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL') &&...
        ST.params.fields.Axisymmetric==0)
    
    RF=ST.params.fields.R;
    %PHIF=ST.params.fields.PHI;
    ZF=ST.params.fields.Z;
    
    %PSIP=zeros(size(RF,1),size(PHIF,1),size(ZF,1));
    
    filename='/Users/21b/Desktop/M3D-C1/it1ip1_p32.rs16/I_psi_0.h5';
    PSI=h5read(filename,'/PSI');
    
    PSIP1=squeeze(PSI(1,:,:));
    FLAG=squeeze(ST.params.fields.Flag(:,1,:));
    
    PHI=mod(PHI,2*pi);
    
    %[RR,PP,ZZ]=meshgrid(RF,PHIF,ZF);
    [RR,ZZ]=meshgrid(RF,ZF);
    %PSIp_part=interp3(RR,PP,ZZ,PSIP1,R,PHI,Z);
    PSIp_part=interp2(RR,ZZ,PSIP1',R,Z);
    
    momgcb=-qe*PSIp_part;
    momgc=momgcv+momgcb;
    
    PSIP=PSIP1;
    
elseif (strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL') &&...
        ST.params.fields.Axisymmetric==1 && isfield(ST.params.fields,'psi_p'))
    
    PSIP=ST.params.fields.psi_p;
    FLAG=ST.params.fields.Flag2D;

    RF=ST.params.fields.R;
    ZF=ST.params.fields.Z;
    
elseif (strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL') &&...
        ST.params.fields.Axisymmetric==1 && isfield(ST.params.fields,'psi_p3D'))
    
    PSIP=ST.params.fields.psi_p3D;
    FLAG=ST.params.fields.Flag3D;

    RF=ST.params.fields.R;
    TF=ST.params.fields.PHI;
    ZF=ST.params.fields.Z; 
    
end


    
%% plot data
%close all

histpeta=0;
histreta=0;
histp=0;
histpmulti=0;
histxi=0;
histRZ=0;
evoEave=0;
evoI=0;
evomomgc=0;
momgccomp=0;
peta_RE_combo=1;
multiRZ=0;
NREevo=0;
fRE=0;

timeslice=size(time,2);

%figure;
%hold on
%scatter(R(:,timeslice),Z(:,timeslice),'o')
%hold off
%daspect([1,1,1])

if fRE==1

    fig=figure;
    
    KEp=me*c^2*(gam-1)/qe;
    
    KEbins=10.^(linspace(log10(min(KEp/(10^6))),log10(max(KEp/(10^6))),50));
    
    [histKE,edgesKE]=histcounts(KEp(flagCon(:,timeslice)>0,timeslice)./(10^6),...
        KEbins);

    hist_binsKE=zeros(size(histKE));
    for i=1:size(histKE,2)
        hist_binsKE(i)=(edgesKE(i)+edgesKE(i+1))/2;
    end

    binwidth=diff(edgesKE);
    normKE=trapz(hist_binsKE,histKE);

    
    loglog(hist_binsKE,histKE/normKE./binwidth,'Linewidth',2)      
    xlabel('E [MeV]')
    ylabel('f_E')
    
    fig=figure;
    loglog(hist_binsKE,histKE,'Linewidth',2)      
    xlabel('E [MeV]')
    ylabel('N_{RE}')
end

if NREevo==1                
    
    E0=ST.params.fields.Eo;
    %E_CH=ST.params.collisions_ss.Ec;
    %tau_c=ST.params.collisions_ss.Tau;
    %Clog0=ST.params.collisions_ss.Clogee;
    %Zeff_c=ST.params.collisions_ss.Zeff;
    %dt_c=ST.params.simulation.dt*ST.params.collisions_ss.subcycling_iterations;
    %ne_c=ST.params.collisions_ss.ne;
    
    
    RPgrowthrate=((E0/E_CH)-1)/(2*tau_c*Clog0);
    
    NREtot=sum(flagRE);
    NREact=sum(flagRE.*flagCol);
    
    startind=int64(.5*size(time,2));
    %startind=1;
    
    myExp='a*exp(b*x)';
    startPointsE=[NREtot(startind),RPgrowthrate];
    startPointsA=[NREact(startind),RPgrowthrate];
      
  
    f1=fit(time(startind:end)',NREtot(startind:end)',myExp,'Start',startPointsE);
    f2=fit(time(startind:end)',NREact(startind:end)',myExp,'Start',startPointsA);
    
    fig=figure;
    p1=semilogy(time,NREtot,'linewidth',3,'color',[0,0,1]);
    hold on
    p1a=semilogy(time,f1.a*exp(time*f1.b),'--','linewidth',3,'color',[0,0,1]);
    p2=semilogy(time,NREact,'linewidth',3,'color',[1,0,0]);
    p2a=semilogy(time,f2.a*exp(time*f2.b),'--','linewidth',3,'color',[1,0,0]);
    
    legend([p1,p1a,p2,p2a],{'All RE','All RE fit','Active RE','Active RE fit'},...
        'numcolumns',2,'Location','southeast')
    
    xlabel('Time(s)')
    ylabel('N_{\rm RE}')
    
    %txt={strcat('Zeff: ',num2str(Zeff_c)),...
    %    strcat('dt=',num2str(dt_c),'(s)'),...
    %    strcat('n_e= ',num2str(ne_c),'m^{-3}'),...
    %    strcat('E/E_{CH}: ',num2str(E0/E_CH)),...        
    %    strcat('\Gamma_{\rm RP}: ',num2str(RPgrowthrate)),...
    %    strcat('All RE fit \gamma: ',num2str(f1.b)),...
    %    strcat('Active RE fit \gamma: ',num2str(f2.b))};
    %text([.025],[.7],txt,'FontSize',14,'Units','normalized')
    
    set(gca,'Fontsize',20)
    
    print(fig,'NREevo','-depsc')
    
end
    
if multiRZ==1
    

    
    fig=figure('Renderer', 'painters', 'Position', [10 10 700 400]);

    
    
    subplot('Position',[.05, .1, .2, .65])
    box on
    hold on    
    if strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL')
        contour(RF,ZF,PSIP',20,'k','LineWidth',1)
    end

    if strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL')
        contour(RF,ZF,FLAG',1,'k','LineWidth',2)
    end
           
    hist0=histogram2(Rj(flag(:,timeslice)>0,timeslice),Zj(flag(:,timeslice)>0,timeslice),30,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    hold off

    set(gca,'Fontsize',12) 

    cb=colorbar();
    cb.Label.String='$N_{\rm RE}$';
    cb.Label.Interpreter='latex';
    cb.Label.FontSize=14;       
    
    xlabel('$R\,{\rm (m)}$','Interpreter','latex','FontSize',14)
    ylabel('$Z\,{\rm (m)}$','Interpreter','latex','FontSize',14)
    title('$\psi_{N,{\rm max}}=0.8446$','interpreter','latex')
    daspect([1,1,1])
    yticks([-1,0,1]);

    colormap(gca,'jet')
    
    subplot('Position',[.275, .1, .2, .65])
    box on
    hold on    
    if strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL')
        contour(RF,ZF,PSIP',20,'k','LineWidth',1)
    end

    if strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL')
        contour(RF,ZF,FLAG',1,'k','LineWidth',2)
    end
           
    hist2=histogram2(Rj2(flag(:,timeslice)>0,timeslice),Zj2(flag(:,timeslice)>0,timeslice),...
        'XBinEdges',hist0.XBinEdges,'YBinEdges',hist0.YBinEdges,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    hold off

    set(gca,'Fontsize',12) 
    

    cb=colorbar();
    cb.Label.String='$N_{\rm RE}$';
    cb.Label.Interpreter='latex';
    cb.Label.FontSize=14;       
    
    xlabel('$R\,{\rm (m)}$','Interpreter','latex','FontSize',14)
    %ylabel('$Z\,{\rm (m)}$','Interpreter','latex','FontSize',14)
    title('$\psi_{N,{\rm max}}=0.63345$','interpreter','latex')
    daspect([1,1,1])
    %yticks([-1,0,1]);
    yticklabels([]);
    
    colormap(gca,'jet')
    
    subplot('Position',[.5, .1, .2, .65])
    box on
    hold on    
    if strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL')
        contour(RF,ZF,PSIP',20,'k','LineWidth',1)
    end

    if strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL')
        contour(RF,ZF,FLAG',1,'k','LineWidth',2)
    end
           
    hist1=histogram2(Rj1(flag(:,timeslice)>0,timeslice),Zj1(flag(:,timeslice)>0,timeslice),...
        'XBinEdges',hist0.XBinEdges,'YBinEdges',hist0.YBinEdges,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    hold off

    set(gca,'Fontsize',12) 
    

    cb=colorbar();
    cb.Label.String='$N_{\rm RE}$';
    cb.Label.Interpreter='latex';
    cb.Label.FontSize=14;       
    
    xlabel('$R\,{\rm (m)}$','Interpreter','latex','FontSize',14)
    %ylabel('$Z\,{\rm (m)}$','Interpreter','latex','FontSize',14)
    title('$\psi_{N,{\rm max}}=0.4223$','interpreter','latex')
    daspect([1,1,1])
    %yticks([-1,0,1]);
    yticklabels([]);
    
    colormap(gca,'jet')
    
    subplot('Position',[.725, .1, .2, .65])
    box on
    hold on    
    if strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL')
        contour(RF,ZF,PSIP',20,'k','LineWidth',1)
    end

    if strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL')
        contour(RF,ZF,FLAG',1,'k','LineWidth',2)
    end
           
    hist3=histogram2(Rj3(flag(:,timeslice)>0,timeslice),Zj3(flag(:,timeslice)>0,timeslice),...
        'XBinEdges',hist0.XBinEdges,'YBinEdges',hist0.YBinEdges,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    hold off

    set(gca,'Fontsize',12) 
    

    cb=colorbar();
    cb.Label.String='$N_{\rm RE}$';
    cb.Label.Interpreter='latex';
    cb.Label.FontSize=14;       
    
    xlabel('$R\,{\rm (m)}$','Interpreter','latex','FontSize',14)
    %ylabel('$Z\,{\rm (m)}$','Interpreter','latex','FontSize',14)
    title('$\psi_{N,{\rm max}}=0.21115$','interpreter','latex')
    daspect([1,1,1])
    %yticks([-1,0,1]);
    yticklabels([]);

    colormap(gca,'jet')
    
    orient(fig,'landscape')
    print(fig,'D3D164409_initspatRZ','-dpdf')
    
end

if histpmulti==1
    fig=figure;

    pmax=20;
    
    timeslices=[1 20 40 60 80];
    pedges=linspace(min(pmag1(flag(:,timeslices(end))>0,...
        timeslices(end))/(me*c)),pmax,100);
    
    hold on
    
    cm=colormap(jet(size(timeslices,2)));
    legends={size(timeslices,2)};
    
    for ii=1:size(timeslices,2) 
        
        [histp,edgesp]=histcounts(pmag(flag(:,timeslices(ii))>0,timeslices(ii))/(me*c)...
            ,pedges,'Normalization','count');

        hist_binsp=zeros(size(histp));
        for i=1:size(histp,2)
            hist_binsp(i)=(edgesp(i)+edgesp(i+1))/2;
        end
        
        legends{ii}=strcat('t=',num2str(time(ii)),'s');
        
        plot(hist_binsp,histp,'Linewidth',2,'color',cm(ii,:))        
    end

    
    set(gca,'Fontsize',12)
    
    %text(-.25,.98,'a)','FontSize',14,'Units','normalized');
    
    ylabel('$N_{\rm RE}$','Interpreter','latex','FontSize',16)
    xlabel('$p$ ($m_ec$)','Interpreter','latex','FontSize',16)
%    axis([0,55,hist.YBinLimits(1),hist.YBinLimits(2)])
    hold off
    
    legend(legends,'location','northwest','fontsize',14)

    print(fig,'TEST25e_histp_evo','-dpdf')
end

if peta_RE_combo==1
    fig=figure('Renderer', 'painters', 'Position', [10 10 750 400]);
    fig.PaperUnits='points';
    fig.PaperPosition = [0 50 650 350];

    
    pfilter=pmag(:,timeslice);
    etafilter=eta(:,timeslice);
    Rfilter=R(:,timeslice);
    Zfilter=Z(:,timeslice);
    
    zmax=0.1;
    zmin=-0.1;
    rmax=1.1;
    rmin=1.0;
    
    Zfilter(Rfilter>rmax)=0;
    Zfilter(Rfilter<rmin)=0;
    Rfilter(Zfilter>zmax)=0;
    Rfilter(Zfilter<zmin)=0;
    Rfilter(Rfilter>rmax)=0;
    Rfilter(Rfilter<rmin)=0;
    Zfilter(Zfilter>zmax)=0;
    Zfilter(Zfilter<zmin)=0;
    
    Rfilter=Rfilter(Rfilter~=0);
    Zfilter=Zfilter(Zfilter~=0);
    pfilter=pfilter(Rfilter~=0);
    etafilter=etafilter(Zfilter~=0);
    
    subplot(1,2,1)   
%    hist=histogram2(etafilter,pfilter/(me*c),50,...
%        'DisplayStyle','tile','LineStyle','none',...
%        'Normalization','count');
    hist=histogram2(pi/180*eta(flagCon(:,timeslice)>0,timeslice),...
        me*c^2*(gam(flagCon(:,timeslice)>0,timeslice)-1)/(10^6*qe),50,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    colormap jet
    cb=colorbar();
    cb.Label.String='$N_{\rm RE}$';
    cb.Label.Interpreter='latex';
    cb.Label.FontSize=14;
    
    set(gca,'Fontsize',12)
    
    text(-.2,.98,'a)','FontSize',14,'Units','normalized');
    
    xlabel('$\eta\,({\rm rad})$','Interpreter','latex','FontSize',16)
    ylabel('$E\,({\rm MeV})$','Interpreter','latex','FontSize',16)
%    axis([0,55,hist.YBinLimits(1),hist.YBinLimits(2)])
    hold off
    
    subplot(1,2,2)
    
    if strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL')
        contour(RF,ZF,squeeze(PSIP(:,11,:))',50,'k','LineWidth',1)
    end
    hold on
    if strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL')
        contour(RF,ZF,squeeze(FLAG(:,11,:))',1,'k','LineWidth',2)
    end
       
    
%    hist=histogram2(Rfilter,Zfilter,30,...
%        'DisplayStyle','tile','LineStyle','none',...
%        'Normalization','count');
    hist=histogram2(R(flagCon(:,timeslice)>0,timeslice),Z(flagCon(:,timeslice)>0,timeslice),30,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    colormap jet
    cb=colorbar();
    cb.Label.String='$N_{\rm RE}$';
    cb.Label.Interpreter='latex';
    cb.Label.FontSize=14;
    
    set(gca,'Fontsize',12)
    
    
    text(-.3,.98,'b)','FontSize',14,'Units','normalized');
    
    
    xlabel('$R\,{\rm (m)}$','Interpreter','latex','FontSize',16)
    ylabel('$Z\,{\rm (m)}$','Interpreter','latex','FontSize',16)
    %axis([rmin,rmax,zmin,zmax])
    daspect([1,1,1])
    hold off

    print(fig,'D3D164409_initial_dist','-djpeg')
end

if histpeta==1
    fig=figure;
    hist=histogram2(eta(:,timeslice),pmag(:,timeslice)/(me*c),50,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    colormap jet
    cb=colorbar();
    cb.Label.String='$N_{\rm RE}$';
    cb.Label.Interpreter='latex';
    cb.Label.FontSize=14;
    
    set(gca,'Fontsize',12)
    
    xlabel('Pitch Angle $\eta$','Interpreter','latex','FontSize',16)
    ylabel('$p$ ($m_ec$)','Interpreter','latex','FontSize',16)
%    axis([0,55,hist.YBinLimits(1),hist.YBinLimits(2)])
    hold off
    
%    txt={strcat('E:',num2str(ST.params.fields.Eo)),...
%         strcat('Impurity:',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('Z0= ',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('dt= ',num2str(ST.params.simulation.dt)),...
%         strcat('ne= ',num2str(ST.params.profiles.neo))};
%    text(.25,.9,txt,'FontSize',14,'Units','normalized');
    
    print(fig,'Initial_peta_distribution','-dpdf')
end

if histreta==1
    fig=figure;
    hist=histogram2(R(:,timeslice),eta(:,timeslice),50,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    colormap jet
    colorbar()
    xlabel('$R (m)$','Interpreter','latex','FontSize',14)
    ylabel('Pitch Angle $\theta$','Interpreter','latex','FontSize',14)
%    axis([0,55,hist.YBinLimits(1),hist.YBinLimits(2)])
    hold off
    
    print(fig,'TEST13_V0','-dpdf')
end

if histp==1
    fig=figure;
    hist=histogram(pmag(:,timeslice)/(me*c),50,'Normalization','count');
    hold on
   
%    line([pmag(1,1)/(me*c),pmag(1,1)/(me*c)],[0,max(hist.Values)],...
%        'LineStyle','--','Linewidth',2.,'Color','red')

    xlabel('$p$ ($m_ec$)','Interpreter','latex','FontSize',14)
    ylabel('$f(p)$','Interpreter','latex','FontSize',14)
%    axis([0,21,0,max(hist.Values)])
    hold off
    set(gca,'Yscale','log')
    
    
    txt={strcat('E:',num2str(ST.params.fields.Eo)),...
%         strcat('Impurity:',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('Z0= ',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('dt= ',num2str(ST.params.simulation.dt)),...
         strcat('ne= ',num2str(ST.params.profiles.neo))};
    text(.25,.9,txt,'FontSize',14,'Units','normalized');
    
    print(fig,'TEST13_V0','-dpdf')
end

if histxi==1
    fig=figure;
    hist=histogram(xi(flag(:,timeslice)>0,timeslice),50,'Normalization','count');
    hold on
    set(gca,'Yscale','log')
   
%    line([xi(1,1),xi(1,1)],[0,max(hist.Values)],...
%        'LineStyle','--','Linewidth',2.,'Color','red')

    xlabel('$\xi$','Interpreter','latex','FontSize',14)
    ylabel('$f(\xi)$','Interpreter','latex','FontSize',14)
%    axis([-1.1,1.1,0,inf])
    hold off
    
    txt={strcat('E:',num2str(ST.params.fields.Eo)),...
%         strcat('Impurity:',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('Z0= ',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('dt= ',num2str(ST.params.simulation.dt)),...
         strcat('ne= ',num2str(ST.params.profiles.neo))};
    text(.25,.9,txt,'FontSize',14,'Units','normalized');
    
    print(fig,'TEST13_V0','-dpdf')
end

if histRZ==1
    fig=figure;
    
    if isfield(ST.params.fields,'psi_p3D')
        t0=1.405;
        tt=t0+time(timeslice);

        [TT,RR,ZZ]=meshgrid(TF,RF,ZF);

        RRF=squeeze(RR(:,1,:));
        ZZF=squeeze(ZZ(:,1,:));
        TTF=tt*ones(size(RRF));

        PSIPtt=interp3(TT,RR,ZZ,PSIP,TTF,RRF,ZZF);
        FLAGtt=interp3(TT,RR,ZZ,FLAG,TTF,RRF,ZZF);
    else
        PSIPtt=PSIP;
        FLAGtt=FLAG;
    end
    
    if strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL')
        contour(RF,ZF,PSIPtt',50,'k','LineWidth',1,...
            'LineStyle','--')
    end
    hold on
    if strcmp(ST.params.simulation.field_model{1}(1:8),'EXTERNAL')
        contour(RF,ZF,FLAGtt',1,'k','LineWidth',2)
    end
    hist=histogram2(R(flagCon(:,timeslice)>0 & flagCol(:,timeslice)>0,timeslice),...
        Z(flagCon(:,timeslice)>0 & flagCol(:,timeslice)>0,timeslice),50,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    colormap jet
    cb=colorbar();
    cb.Label.String='$N_{\rm RE}$';
    cb.Label.Interpreter='latex';
    cb.Label.FontSize=14;
    
    set(gca,'Fontsize',12)
    
    xlabel('$R\,{\rm (m)}$','Interpreter','latex','FontSize',16)
    ylabel('$Z\,{\rm (m)}$','Interpreter','latex','FontSize',16)
%    axis([1.4,1.85,-.175,.175])
    daspect([1,1,1])

    print(fig,'Initial_RZ_distribution','-dpdf')
end

if evoEave==1
    fig=figure; 
    plot(time,AveE)
    xlabel('${\rm t\,(s)}$','Interpreter','latex','FontSize',14)
    ylabel('${\rm E_{\rm ave}\,(eV)}$','Interpreter','latex','FontSize',14)

    txt={strcat('E:',num2str(ST.params.fields.Eo)),...
%         strcat('Impurity:',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('Z0= ',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('dt= ',num2str(ST.params.simulation.dt)),...
         strcat('ne= ',num2str(ST.params.profiles.neo))};
    text(.25,.9,txt,'FontSize',14,'Units','normalized');
    
    print(fig,'TEST13_evoEave','-dpdf')
    
end
    
if evoI==1
    fig=figure; 
    plot(time,I/I(1))
    xlabel('${\rm t\,(s)}$','Interpreter','latex','FontSize',14)
    ylabel('${\rm I/I(0)}$','Interpreter','latex','FontSize',14)

    txt={strcat('E:',num2str(ST.params.fields.Eo)),...
%         strcat('Impurity:',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('Z0= ',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('dt= ',num2str(ST.params.simulation.dt)),...
         strcat('ne= ',num2str(ST.params.profiles.neo))};
    text(.25,.9,txt,'FontSize',14,'Units','normalized');
    
    print(fig,'evoI','-dpdf')

end

if evomomgc==1
    fig=figure;

    momnorm=zeros(size(momgc));
    for i=1:size(momgc,1)
        momnorm(i,:)=(momgc(i,:)-momgc(i,1))/momgc(i,1);
    end
    
    avemomnorm=sum(momnorm.*flag)./sum(flag);
    
    plot(time,momnorm(1,:))
    hold on
    line([0,time(end)],[0,0],'Color','k','LineStyle',':')    
    xlabel('${\rm t\,(s)}$','Interpreter','latex','FontSize',14)
    ylabel('${\rm \Delta p_\zeta/p_\zeta(0)}$','Interpreter','latex','FontSize',14)
    hold off

    txt={strcat('dt= ',num2str(ST.params.simulation.dt))};
    text(.25,.9,txt,'FontSize',14,'Units','normalized');
    
    print(fig,'evoI','-dpdf')

end

if momgccomp==1
    fig=figure;


    momnorm0=(momgc0(1,:)-momgc0(1,1))/momgc0(1,1);
    momnorm1=(momgc1(1,:)-momgc1(1,1))/momgc1(1,1);
    momnorm2=(momgc2(1,:)-momgc2(1,1))/momgc2(1,1);
    momnorm3=(momgc3(1,:)-momgc3(1,1))/momgc3(1,1);
    momnorm4=(momgc4(1,:)-momgc4(1,1))/momgc4(1,1);
    momnorm5=(momgc5(1,:)-momgc5(1,1))/momgc5(1,1);
    momnorm6=(momgc6(1,:)-momgc6(1,1))/momgc6(1,1);
    momnorm7=(momgc7(1,:)-momgc7(1,1))/momgc7(1,1);   
    
    semilogy(time0,abs(momnorm0))
    hold on
    semilogy(time1,abs(momnorm1))
    semilogy(time2,abs(momnorm2))
    semilogy(time3,abs(momnorm3))
    semilogy(time4,abs(momnorm4))
    semilogy(time5,abs(momnorm5))
    semilogy(time6,abs(momnorm6))
    semilogy(time7,abs(momnorm7))
   
    xlabel('${\rm t\,(s)}$','Interpreter','latex','FontSize',14)
    ylabel('${\rm p_\zeta/p_\zeta(0)}$','Interpreter','latex','FontSize',14)
    hold off
    legend({'0','1','2','3','4','5','6','7'},'Location','Northwest')

    
    print(fig,'compare_momgct','-dpdf')

end
