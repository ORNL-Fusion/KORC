function ST = diagnoseKORC(path,visible,range)
% ST = diagnoseKORC('../KORC-FO/outputFiles/','on',[0,100])
%close all

ST = struct;
ST.path = path;
ST.visible = visible;

ST.range = range;

% ST.num_snapshots = ST.range(2) - ST.range(1) + 1;

ST.params = loadSimulationParameters(ST);

% ST.time = ...
%     ST.params.simulation.dt*double(ST.params.simulation.output_cadence)*...
%     double(ST.range(1):1:ST.range(2));

ST = loadData(ST);

% energyConservation(ST);

% angularMomentum(ST);

% ST.RT = radialTransport(ST);

% ST.CP = confined_particles(ST);
% 
% ST.PAD = pitchAngleDiagnostic(ST,30);

% ST.MMD = magneticMomentDiagnostic(ST,70);

% poloidalPlaneDistributions(ST,25);

% ST.CMF = changeOfMagneticField(ST);

% energyLimit(ST);

% ST.PR = LarmorVsLL(ST);

% stackedPlots(ST,40);

% scatterPlots(ST);

% ST.P = synchrotronSpectrum(ST,true);

% ST.SD = syntheticDiagnosticSynchrotron(ST,false);

% radiationPlane(ST);

% calculateTemperatureComponents(ST);

% SE_phaseSpaceAnalisys(ST);

% plotEnergyPitchanglePDF(ST);

% figuresAPS2017(ST);

% NIMROD_figure(ST);

% ST.M = movieEnergyPitchAngle(ST);

% pitchAnglePDFSlices(ST)
% 
% PAE_PDF(ST);

% pitchAnglePDFSlicesVsTime(ST)

% spatialAndVelocityPDF(ST)

% subPhaseSpaceAnalysis(ST);

% plotEnergyPDF(ST);

% save('energy_limit','ST')
end

function params = loadSimulationParameters(ST)
params = struct;

info = h5info([ST.path 'simulation_parameters.h5']);

for ii=1:length(info.Groups)
    for jj=1:length(info.Groups(ii).Datasets)
        name = info.Groups(ii).Name(2:end);
        subname = info.Groups(ii).Datasets(jj).Name;
        params.(name).(subname) = ...
            h5read(info.Filename,['/' name '/' subname]);
    end
end

try
    info = h5info([ST.path 'avalanche_parameters.h5']);
    
    for ii=1:length(info.Groups)
        for jj=1:length(info.Groups(ii).Datasets)
            name = info.Groups(ii).Name(2:end);
            subname = info.Groups(ii).Datasets(jj).Name;
            params.(name).(subname) = ...
                h5read(info.Filename,['/' name '/' subname]);
        end
    end
catch
end

try
    info = h5info([ST.path 'simple_equilibrium_pdf.h5']);
    
    for ii=1:length(info.Groups)
        for jj=1:length(info.Groups(ii).Datasets)
            name = info.Groups(ii).Name(2:end);
            subname = info.Groups(ii).Datasets(jj).Name;
            params.(name).(subname) = ...
                h5read(info.Filename,['/' name '/' subname]);
        end
    end
catch
end

ST.range(2)+1;
try
    info = h5info([ST.path 'experimental_distribution_parameters.h5']);
    
    for ii=1:length(info.Groups)
        for jj=1:length(info.Groups(ii).Datasets)
            name = info.Groups(ii).Name(2:end);
            subname = info.Groups(ii).Datasets(jj).Name;
            params.(name).(subname) = ...
                h5read(info.Filename,['/' name '/' subname]);
        end
    end
catch
end

try
    info = h5info([ST.path 'synthetic_camera.h5']);
    for ii=1:length(info.Groups)
        for jj=1:length(info.Groups(ii).Datasets)
            name = info.Groups(ii).Name(2:end);
            subname = info.Groups(ii).Datasets(jj).Name;
            params.(name).(subname) = ...
                h5read(info.Filename,['/' name '/' subname]);
        end
    end
catch
end

end

function ST = loadData(ST)
data = struct;

list = ST.params.simulation.outputs_list;

filename = [ST.path 'file_0.h5'];
H = h5info(filename);
it = zeros(1,size(H.Groups,1));
for ii=1:size(H.Groups,1)
    tmpstr = strrep(H.Groups(ii).Name,'/','');
    it(ii) = str2num(tmpstr);
end


it = sort(it,'ascend');


cadence = double(ST.params.simulation.output_cadence);
[~,I1] = min(abs(cadence*ST.range(1)-it));
[~,I2] = min(abs(cadence*ST.range(2)-it));

indices = it(I1:I2);

ST.num_snapshots = I2 - I1 + 1;

% ST.params = loadSimulationParameters(ST);

ST.time = zeros(1,ST.num_snapshots);

for ll=1:length(list)
    for ss=1:ST.params.simulation.num_species
        tnp = double(ST.params.species.ppp(ss)*ST.params.simulation.nmpi);

        
        if (strcmp(list{ll},'X') || strcmp(list{ll},'Y') || strcmp(list{ll},'V') || ...
                strcmp(list{ll},'B') || strcmp(list{ll},'gradB') || ...
                strcmp(list{ll},'curlb') || strcmp(list{ll},'E') || ...
                strcmp(list{ll},'Y0') || strcmp(list{ll},'Y1') || ...
                strcmp(list{ll},'Yborn'))
            data.(['sp' num2str(ss)]).(list{ll}) = zeros(tnp,3,ST.num_snapshots);
        elseif (strcmp(list{ll},'RHS'))
            data.(['sp' num2str(ss)]).(list{ll}) = zeros(tnp,5,ST.num_snapshots);
        elseif (strcmp(list{ll},'E_SC') || strcmp(list{ll},'J_SC'))
            tnp_SC=size(ST.params.fields.PSIP_1D,1);
            data.(['sp' num2str(ss)]).(list{ll}) = zeros(tnp_SC,ST.num_snapshots);
        else
            data.(['sp' num2str(ss)]).(list{ll}) = zeros(tnp,ST.num_snapshots);            
        end
    end
end

for impi=1:ST.params.simulation.nmpi
    disp(['Loading file:' num2str(impi) ' out of ' num2str(ST.params.simulation.nmpi)])
    filename = [ST.path 'file_' num2str(impi-1) '.h5'];
    for ss=1:ST.params.simulation.num_species
        indi = (impi - 1)*double(ST.params.species.ppp(ss)) + 1;
        indf = impi*double(ST.params.species.ppp(ss));
        
        for ii=1:ST.num_snapshots
            if mod(ii,floor(ST.num_snapshots/10))==0
               disp(ii) 
            end
            for ll=1:length(list)
                try
                    dataset = ...
                        ['/' num2str(indices(ii)) '/spp_' num2str(ss)...
                        '/' list{ll}];

                    if (strcmp(list{ll},'X') || strcmp(list{ll},'Y') ||strcmp(list{ll},'V') || ...
                            strcmp(list{ll},'B') || strcmp(list{ll},'gradB') || ...
                            strcmp(list{ll},'curlb') ||strcmp(list{ll},'E')|| ...
                            strcmp(list{ll},'Y0') || strcmp(list{ll},'Y1') || ...
                            strcmp(list{ll},'Yborn'))
                        
                        test=h5read(filename, dataset);
                        
                        data.(['sp' num2str(ss)]).(list{ll})(indi:indf,:,ii) = test;
                        
                    elseif (strcmp(list{ll},'RHS'))
                        data.(['sp' num2str(ss)]).(list{ll})(indi:indf,:,ii) = ...
                            h5read(filename, dataset);
                    elseif (strcmp(list{ll},'E_SC') || strcmp(list{ll},'J_SC'))
                        
                        test=h5read(filename, dataset);
                        
                        data.(['sp' num2str(ss)]).(list{ll})(:,ii) = test;
                    else
                        data.(['sp' num2str(ss)]).(list{ll})(indi:indf,ii) = ...
                            h5read(filename, dataset);
                    end
                
                    

                    dataset = ...
                        ['/' num2str(indices(ii)) '/time'];
            
                    
                    ST.time(ii) = h5read(filename, dataset);
                    

                    
                catch EX
                    if ~isempty(strfind(EX.message,'not found'))
                        ST.num_snapshots = ST.num_snapshots - 1;
                        ST.time(ii) = [];
                        indices(ii) = [];
                        for jj=1:length(list)
                            if (strcmp(list{jj},'X') || strcmp(list{jj},'Y') || strcmp(list{jj},'V') || strcmp(list{jj},'B') || strcmp(list{ll},'gradB') || strcmp(list{ll},'curlb') ||strcmp(list{jj},'E'))
                                data.(['sp' num2str(ss)]).(list{jj})(:,:,ii) = [];
                            else
                                data.(['sp' num2str(ss)]).(list{jj})(:,ii) = [];
                            end
                        end
                        break
                    end
                end
            end
        end
    end
end


% for ll=1:length(list)
%     disp(['Loading ' list{ll}])
%     for ss=1:ST.params.simulation.num_species
%         for impi=1:ST.params.simulation.nmpi
%             filename = [ST.path 'file_' num2str(impi-1) '.h5'];
%             indi = (impi - 1)*double(ST.params.species.ppp(ss)) + 1;
%             indf = impi*double(ST.params.species.ppp(ss));
%             for ii=I1:I2
%                 try
%                     dataset = ...
%                         ['/' num2str(it(ii)) '/spp_' num2str(ss)...
%                         '/' list{ll}];
%                     if (strcmp(list{ll},'X') || strcmp(list{ll},'V') || strcmp(list{ll},'B') || strcmp(list{ll},'E'))
%                         data.(['sp' num2str(ss)]).(list{ll})(:,indi:indf,ii) = ...
%                             h5read(filename, dataset);
%                     else
%                         data.(['sp' num2str(ss)]).(list{ll})(indi:indf,ii) = ...
%                             h5read(filename, dataset);
%                     end
%                 catch EX
%                     if ~isempty(strfind(EX.message,'not found'))
%                         ST.num_snapshots = ST.num_snapshots - 1;
%                         ST.time(ii) = [];
%                         it(ii) = [];
%                         for jj=1:ST.params.simulation.num_species
%                             if (strcmp(list{jj},'X') || strcmp(list{jj},'V') || strcmp(list{jj},'B') || strcmp(list{jj},'E'))
%                                 data.(['sp' num2str(ss)]).(list{jj})(:,:,ii) = [];
%                             else
%                                 data.(['sp' num2str(ss)]).(list{jj})(:,ii) = [];
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

ST.data = data;

end

function energyConservation(ST)
err = zeros(ST.num_snapshots,ST.params.simulation.num_species);
maxerr = zeros(ST.num_snapshots,ST.params.simulation.num_species);
minerr = zeros(ST.num_snapshots,ST.params.simulation.num_species);

st1 = zeros(ST.num_snapshots,ST.params.simulation.num_species);
st2 = zeros(ST.num_snapshots,ST.params.simulation.num_species);
st3 = zeros(ST.num_snapshots,ST.params.simulation.num_species);
st4 = zeros(ST.num_snapshots,ST.params.simulation.num_species);

h1 = figure('Visible',ST.visible);
h2 = figure('Visible',ST.visible);
h3 = figure('Visible',ST.visible);
h4 = figure('Visible',ST.visible);
h5 = figure('Visible',ST.visible);
h6 = figure('Visible',ST.visible);

set(h1,'name','Energy conservation','numbertitle','off')
set(h2,'name','Velocity components','numbertitle','off')
set(h3,'name','Radiated power','numbertitle','off')
set(h4,'name','Input power','numbertitle','off')
set(h5,'name','Energy gain/loss','numbertitle','off')
set(h6,'name','Energy statistics','numbertitle','off')
% try
    for ss=1:ST.params.simulation.num_species
        m = ST.params.species.m(ss);
        q = abs(ST.params.species.q(ss));
        c = ST.params.scales.v;
        
        
        pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
%         pin = true(1,size(ST.data.(['sp' num2str(ss)]).g,1));
%         passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
%         bool = pin & passing;

        gammap = ST.data.(['sp' num2str(ss)]).g(pin,:);
        tmp = zeros(size(gammap));
        for ii=1:size(tmp,1)
            tmp(ii,:) = 100*( gammap(ii,:) - gammap(ii,1) )./gammap(ii,1);
%             tmp(ii,:) = (gammap(ii,:)*m*c^2/q)/1E6;
        end
        
        err(:,ss) = mean(tmp,1);
        maxerr(:,ss) = max(tmp,[],1);
        minerr(:,ss) = min(tmp,[],1);
%         maxerr(:,ss) = err(:,ss) + std(tmp,0,1)';
%         minerr(:,ss) = err(:,ss) - std(tmp,0,1)';
        
        if (~isempty(tmp))
            st1(:,ss) = mean(tmp,1);
            st2(:,ss) = std(tmp,0,1);
            st3(:,ss) = zeros(size(st1(:,ss)));%skewness(tmp,1,1);
            st4(:,ss) = zeros(size(st1(:,ss)));%kurtosis(tmp,1,1);
        end
        
        try
            Prad = mean(abs(ST.data.(['sp' num2str(ss)]).Prad(pin,:)),1);
            minPrad = Prad + std(abs(ST.data.(['sp' num2str(ss)]).Prad(pin,:)),0,1);
            maxPrad = Prad - std(abs(ST.data.(['sp' num2str(ss)]).Prad(pin,:)),0,1);
        catch
            Prad = zeros(size(ST.time));
            minPrad = 0;
            maxPrad = 1;
        end
        
        
        try
            Pin = mean(abs(ST.data.(['sp' num2str(ss)]).Pin(pin,:)),1);
            minPin = Pin + std(abs(ST.data.(['sp' num2str(ss)]).Pin(pin,:)),0,1);
            maxPin = Pin - std(abs(ST.data.(['sp' num2str(ss)]).Pin(pin,:)),0,1);
        catch
            Pin = zeros(size(ST.time));
            minPin = 0;
            maxPin = 1;
        end
        
        
        eta = deg2rad(ST.data.(['sp' num2str(ss)]).eta(pin,:));
        V = sqrt(1 - 1./gammap.^2); % in units of the speed of light
        Vpar = mean(V.*cos(eta),1);
        Vperp = mean(V.*sin(eta),1);
        
        figure(h1)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        try
            plot(ST.time,err(:,ss),'k-',ST.time,minerr(:,ss),'r:',ST.time,maxerr(:,ss),'r:')
%             hold on
%             plot(ST.time,tmp)
%             hold off
        catch
        end
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\Delta \mathcal{E}/\mathcal{E}_0$ (\%)','Interpreter','latex','FontSize',16)
        saveas(h1,[ST.path 'energy_conservation'],'fig')
        
        figure(h2)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(ST.time,Vpar./Vperp,'k-')
%         plot(time,Vpar,'k-',time,Vperp,'r-')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$v_\parallel/v_\perp$ ($c$)','Interpreter','latex','FontSize',16)
        saveas(h2,[ST.path 'velocity_components'],'fig')

        figure(h6)
        subplot(4,1,1)
        hold on; plot(ST.time,st1(:,ss)); hold off
        box on; grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\mu$','Interpreter','latex','FontSize',16)
        subplot(4,1,2)
        hold on; plot(ST.time,st2(:,ss)); hold off
        box on; grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\sigma$','Interpreter','latex','FontSize',16)
        subplot(4,1,3)
        hold on; plot(ST.time,st3(:,ss)); hold off
        box on; grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$s$','Interpreter','latex','FontSize',16)
        subplot(4,1,4)
        hold on; plot(ST.time,st4(:,ss)); hold off
        box on; grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$k$','Interpreter','latex','FontSize',16)
        saveas(h6,[ST.path 'energy_statistics'],'fig')
        
        
        figure(h3)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(ST.time,Prad,'k',ST.time,minPrad,'r:',ST.time,maxPrad,'r:')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\langle P_{rad} \rangle$','Interpreter','latex','FontSize',16)
        saveas(h3,[ST.path 'radiated_power'],'fig')
        
        figure(h4)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(ST.time,Pin,'k',ST.time,minPin,'r:',ST.time,maxPin,'r:')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\langle P_{in} \rangle$','Interpreter','latex','FontSize',16)
        saveas(h4,[ST.path 'input_power'],'fig')
        
        figure(h5)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(ST.time,Prad./Pin,'k')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\langle P_{rad} \rangle / \langle P_{in} \rangle$','Interpreter','latex','FontSize',16)
        saveas(h5,[ST.path 'energy_gain-lost'],'fig')
    end
% catch
%     error('Something went wrong: energyConservation')
% end
end

function PAD = pitchAngleDiagnostic(ST,nbins)
PAD = struct;
N = 10;

tmax = max(ST.time);
tmin = min(ST.time);

mean_fx = zeros(ST.params.simulation.num_species,ST.num_snapshots);
std_fx = zeros(ST.params.simulation.num_species,ST.num_snapshots);
skewness_fx = zeros(ST.params.simulation.num_species,ST.num_snapshots);
kurtosis_fx = zeros(ST.params.simulation.num_species,ST.num_snapshots);

mean_fz = zeros(ST.params.simulation.num_species,ST.num_snapshots);
std_fz = zeros(ST.params.simulation.num_species,ST.num_snapshots);
skewness_fz = zeros(ST.params.simulation.num_species,ST.num_snapshots);
kurtosis_fz = zeros(ST.params.simulation.num_species,ST.num_snapshots);

stats = zeros(2,ST.params.simulation.num_species,ST.num_snapshots);

fx = zeros(ST.params.simulation.num_species,nbins,ST.num_snapshots);
x = zeros(ST.params.simulation.num_species,nbins,ST.num_snapshots);
fz = zeros(ST.params.simulation.num_species,nbins,ST.num_snapshots);
z = zeros(ST.params.simulation.num_species,nbins,ST.num_snapshots);

if ST.visible
    h1 = figure('Visible','on');
else
    h1 = figure('Visible','off');
end

set(h1,'name','Pitch angle PDF vs. time','numbertitle','off')
% h0 = figure('Visible',ST.visible);
% set(h0,'name','Energy PDF vs. time','numbertitle','off')
if ST.visible
    h = figure('Visible','on');
else
    h = figure('Visible','off');
end
set(h,'name','Pitch angle: variability','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    q = abs(ST.params.species.q(ss));
    m = ST.params.species.m(ss);
    c = ST.params.scales.v;
    
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    bool = pin; % & passing;
    eta = ST.data.(['sp' num2str(ss)]).eta(bool,:);
    Eo = ST.data.(['sp' num2str(ss)]).g(bool,:)*m*c^2/q;
    Eo = Eo/1E6;
    
%     try
%         tmp = eta(:,end);
%         save(['pitch_' num2str(ST.params.species.etao(ss)) '_E_' num2str(floor(Eo(1,1)))],'tmp')
%     catch
%     end
    
    if ~isempty(eta)
        mean_fx(ss,:) = mean(eta,1);
        std_fx(ss,:) = std(eta,0,1);
        skewness_fx(ss,:) = zeros(size(mean_fx(ss,:)));%skewness(eta,0,1);
        kurtosis_fx(ss,:) = zeros(size(mean_fx(ss,:)));%kurtosis(eta,1,1);
        
        mean_fz(ss,:) = mean(Eo,1);
        std_fz(ss,:) = std(Eo,0,1);
        skewness_fz(ss,:) = zeros(size(mean_fz(ss,:)));%skewness(Eo,0,1);
        kurtosis_fz(ss,:) = zeros(size(mean_fz(ss,:)));%kurtosis(Eo,1,1);
        
        for ii=1:ST.num_snapshots
            minVal = min(min( eta(:,ii) ));
            maxVal = max(max( eta(:,ii) ));
            x(ss,:,ii) = linspace(minVal,maxVal,nbins);
            Dx = mean(diff(x(ss,:,ii)));
            
            minVal = min(min( Eo(:,ii) ));
            maxVal = max(max( Eo(:,ii) ));
            z(ss,:,ii) = linspace(minVal,maxVal,nbins);
            Dz = mean(diff(z(ss,:,ii)));
            
            stats(1,ss,ii) = 100*mean(abs(eta(:,ii) - eta(:,1))./eta(:,1));
            stats(2,ss,ii) = 100*std(abs(eta(:,ii) - eta(:,1))./eta(:,1));
            [fx(ss,:,ii),~] = hist(eta(:,ii),x(ss,:,ii));
            fx(ss,:,ii) = fx(ss,:,ii)/(Dx*sum(fx(ss,:,ii)));
            [fz(ss,:,ii),~] = hist(Eo(:,ii),z(ss,:,ii));
            fz(ss,:,ii) = fz(ss,:,ii)/(Dz*sum(fz(ss,:,ii)));
        end
    end
    
    figure(h1)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    surf(ST.time,squeeze(x(ss,:,:)),log10(squeeze(fx(ss,:,:))),'LineStyle','none')
%     surf(ST.time,squeeze(x(ss,:)),squeeze(fx(ss,:,:)),'LineStyle','none')
%     axis([tmin tmax minVal maxVal])
    box on
    axis on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$\theta$ ($^\circ$)','Interpreter','latex','FontSize',16)
    colormap(jet(256))
    
    figure(h)
    subplot(2,1,1)
    hold on
    plot(ST.time,squeeze(stats(1,ss,:)))
    hold off
    xlim([tmin tmax])
    box on
    axis on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('mean($\Delta |\theta|/\theta_0$) ($\%$)','Interpreter','latex','FontSize',16)
    subplot(2,1,2)
    hold on
    plot(ST.time,squeeze(stats(2,ss,:)))
    hold off
    xlim([tmin tmax])
    box on
    axis on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('std($\Delta |\theta|/\theta_0$) ($\%$)','Interpreter','latex','FontSize',16)
end

if ST.visible
    h2 = figure('Visible','on');
else
    h2 = figure('Visible','off');
end
set(h2,'name','Statistical moments pitch angle','numbertitle','off')
for ii=1:ST.params.simulation.num_species
    figure(h2)
    subplot(4,1,1)
    hold on
    plot(ST.time,mean_fx(ii,:))
    hold off
    figure(h2)
    subplot(4,1,2)
    hold on
    plot(ST.time,std_fx(ii,:))
    hold off
    figure(h2)
    subplot(4,1,3)
    hold on
    plot(ST.time,skewness_fx(ii,:))
    hold off
    figure(h2)
    subplot(4,1,4)
    hold on
    plot(ST.time,kurtosis_fx(ii,:))
    hold off
end

figure(h2)
subplot(4,1,1)
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('mean($\theta$)','Interpreter','latex','FontSize',16)
box on
grid on
figure(h2)
subplot(4,1,2)
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('std($\theta$)','Interpreter','latex','FontSize',16)
box on
grid on
figure(h2)
subplot(4,1,3)
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('skewness($\theta$)','Interpreter','latex','FontSize',16)
box on
grid on
figure(h2)
subplot(4,1,4)
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('kurtosis($\theta$)','Interpreter','latex','FontSize',16)
box on
grid on

offset = floor(double(ST.num_snapshots)/N);

% z = linspace(-4,4,100);
% fz = exp( -0.5*z.^2 )/sqrt(2*pi);

if ST.visible
    h3 = figure('Visible','on');
else
    h3 = figure('Visible','off');
end
set(h3,'name','PDF pitch angle','numbertitle','off')
if ST.visible
    hh = figure('Visible','on');
else
    hh = figure('Visible','off');
end
set(hh,'name','PDF energy','NumberTitle','off')

for ii=1:N
    it = ii*offset;
    nc = floor(N/2);
    nr = floor(N/nc);
    figure(h3)
    subplot(nr,nc,ii)
    hold on
    figure(hh)
    subplot(nr,nc,ii)
    hold on
    
    for ss=1:ST.params.simulation.num_species    
%         xAxis = ( x(ss,:,it) - mean(x(ss,:,it)) )/std(fx(ss,:,it));
%         f = std(fx(ss,:,it))*squeeze(fx(ss,:,it))/trapz(x(ss,:,it),fx(ss,:,it));
        xAxis = x(ss,:,it);
        f = squeeze(fx(ss,:,it));
        
        figure(h3)
        plot(xAxis,f,'o:')
%         plot(xAxis,log10(f),'o:')
                
        xAxis = z(ss,:,it);
        f = squeeze(fz(ss,:,it));
        
        figure(hh)
        plot(xAxis,f,'o:')
    end
    figure(h3)
    title(['Time: ' num2str(ST.time(it))],'Interpreter','latex','FontSize',11)
    hold off
    box on
    grid on
    figure(hh)
    title(['Time: ' num2str(ST.time(it))],'Interpreter','latex','FontSize',11)
    hold off
    box on
    grid on
end

saveas(h,[ST.path 'variability'],'fig')
% saveas(h0,[ST.path 'Energy_PDF_vs_time'],'fig')
saveas(h1,[ST.path 'pitch_vs_time'],'fig')
saveas(h2,[ST.path 'pitch_stats'],'fig')
saveas(h3,[ST.path 'pitch_pdfs'],'fig')
saveas(hh,[ST.path 'energy_pdfs'],'fig')

PAD.mean = mean_fx;
PAD.std = std_fx;

end

function MMD = magneticMomentDiagnostic(ST,numBins)
MMD = struct;
tmax = max(ST.time);
tmin = min(ST.time);

mean_f = zeros(ST.params.simulation.num_species,ST.num_snapshots);
std_f = zeros(ST.params.simulation.num_species,ST.num_snapshots);
skewness_f = zeros(ST.params.simulation.num_species,ST.num_snapshots);
kurtosis_f = zeros(ST.params.simulation.num_species,ST.num_snapshots);

fx = zeros(ST.params.simulation.num_species,numBins,ST.num_snapshots);
x = zeros(ST.params.simulation.num_species,numBins);

if ST.visible
    h1 = figure('Visible','on');
else
    h1 = figure('Visible','off');
end
set(h1,'name','PDF of magnetic moment','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    tmp = ST.data.(['sp' num2str(ss)]).mu(pin,:);
    
    if ~isempty(tmp)
        for ii=1:size(tmp,1)
            tmp(ii,:) = 100*(tmp(ii,:) - tmp(ii,1))./tmp(ii,1);
        end
        
        mean_f(ss,:) = mean(tmp,1);
        std_f(ss,:) = std(tmp,0,1);
        skewness_f(ss,:) = skewness(tmp,0,1);
        kurtosis_f(ss,:) = kurtosis(tmp,1,1);
        
        minVal = min(min( tmp ));
        maxVal = max(max( tmp ));
        x(ss,:) = linspace(minVal,maxVal,numBins);
        
        for ii=1:ST.num_snapshots
            try
                [fx(ss,:,ii),~] = hist(tmp(:,ii),x(ss,:));
                fx(ss,:,ii) = fx(ss,:,ii)/( trapz(x(ss,:),fx(ss,:,ii)) );
%                 fx(ss,:,ii) = fx(ss,:,ii)/max(squeeze(fx(ss,:,ii)));
            catch
            end
        end
    end
    
    subplot(double(ST.params.simulation.num_species),1,double(ss))
%     surf(ST.time,squeeze(x(ss,:)),log10(squeeze(fx(ss,:,2:end))),'LineStyle','none')
    surf(ST.time(2:end),squeeze(x(ss,:)),squeeze(fx(ss,:,2:end)),'LineStyle','none')
%     axis([tmin tmax minVal maxVal])
    view([0,90])
    box on
    axis on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$\Delta \mu/\mu_0$ ($\%$)','Interpreter','latex','FontSize',16)
    colormap(jet(1024))
end

if ST.visible
    h2 = figure('Visible','on');
else
    h2 = figure('Visible','off');
end
set(h2,'name','Statistical moments magnetic moment','numbertitle','off')
for ii=1:ST.params.simulation.num_species
    figure(h2)
    subplot(4,1,1)
    hold on
    plot(ST.time,mean_f(ii,:))
    hold off
    figure(h2)
    subplot(4,1,2)
    hold on
    plot(ST.time,std_f(ii,:))
    hold off
    figure(h2)
    subplot(4,1,3)
    hold on
    plot(ST.time,skewness_f(ii,:))
    hold off
    figure(h2)
    subplot(4,1,4)
    hold on
    plot(ST.time,kurtosis_f(ii,:))
    hold off
end

saveas(h1,[ST.path 'magnetic_moment_pdfs'],'fig')
saveas(h2,[ST.path 'magnetic_moment_stats'],'fig')

MMD.fx = fx;
MMD.x = x;
MMD.stat1 = mean_f;
MMD.stat2 = std_f;
MMD.stat3 = skewness_f;
MMD.stat4 = kurtosis_f;
end

function poloidalPlaneDistributions(ST,nbins)
N = 3;
offset = floor(double(ST.num_snapshots)/N);

for ss=1:ST.params.simulation.num_species
    for ii=1:N
        it = ii*offset;
        
        R = squeeze( sqrt( ST.data.(['sp' num2str(ss)]).X(1,:,it).^2 + ...
            ST.data.(['sp' num2str(ss)]).X(2,:,it).^2 ) );
        Z = squeeze( ST.data.(['sp' num2str(ss)]).X(3,:,it) );
        Prad = squeeze( ST.data.(['sp' num2str(ss)]).Prad(:,it) );

        % Poloidal distribution of particles
        h = figure;
        subplot(3,1,1)
        n = histogram2(R,Z,[nbins,nbins],'FaceColor','flat','Normalization','probability');
        colorbar
        xAxis = n.XBinEdges;
        dx = mean( diff(xAxis) );
        yAxis = n.YBinEdges;
        dy = mean( diff(yAxis) );
        axis([min(xAxis) max(xAxis) min(yAxis) max(yAxis)])
        axis equal
        view([0 90])
        colormap(jet(1024))
        title(['Species: ' num2str(ss) 'Time: ' num2str(ST.time(it))],'Interpreter','latex','FontSize',11)
        xlabel('$R$','Interpreter','latex','FontSize',16)
        ylabel('$Z$','Interpreter','latex','FontSize',16)
        
        x = ST.params.fields.Ro + ...
            ST.params.species.ro(ss)*cos(linspace(0,2*pi,100));
        y = ST.params.species.ro(ss)*sin(linspace(0,2*pi,100));
        hold on; plot3(x,y,1*ones(size(x)),'k');hold off
        
        
        % Poloidal distribution of radiated synchroton power
%        prad = histogram2(R,Z,[nbins,nbins],'FaceColor','flat')       
        prad = zeros(nbins,nbins);
        
        R = R - min(xAxis);
        Z = Z + abs(min(yAxis));
        indx = floor( (R + 0.5*dx)/dx ) + 1;
        indy = floor( (Z + 0.5*dy)/dy ) + 1;
        
        I = find(indx > nbins);
        indx(I) = indx(I) - 1;
        
        I = find(indy > nbins);
        indy(I) = indy(I) - 1;
        
        for jj=1:ST.params.species.ppp(ss)
            prad(indy(jj),indx(jj)) = prad(indx(jj),indy(jj)) + Prad(jj);
        end
        
        figure(h)
        subplot(3,1,2)
        surf(xAxis(1:end-1)',yAxis(1:end-1)',prad,'LineStyle','none')
        axis equal
        view([0 90])        
        colormap(jet(1024))
        axis([min(xAxis) max(xAxis) min(yAxis) max(yAxis)])
        xlabel('$R$','Interpreter','latex','FontSize',16)
        ylabel('$Z$','Interpreter','latex','FontSize',16)
        
        A = n.Values'.*prad;
%         I = isinf(A);
%         A(I) = 0;
        
        figure(h)
        subplot(3,1,3)
        surf(xAxis(1:end-1)',yAxis(1:end-1)',A,'LineStyle','none')
        axis equal
        view([0 90])        
        colormap(jet(1024))
        axis([min(xAxis) max(xAxis) min(yAxis) max(yAxis)])
        xlabel('$R$','Interpreter','latex','FontSize',16)
        ylabel('$Z$','Interpreter','latex','FontSize',16)

    end

% for ii=1:N
%     it = ii*offset;
%     
%     R = squeeze( sqrt( ST.data.(['sp' num2str(ss)]).X(1,:,it).^2 + ...
%         ST.data.(['sp' num2str(ss)]).X(2,:,it).^2 ) );
%     Z = squeeze( ST.data.(['sp' num2str(ss)]).X(3,:,it) );
%     Prad = squeeze( ST.data.(['sp' num2str(ss)]).Prad(:,it) );
%     
%     % Poloidal distribution of particles
%     figure(m)
%     n = histogram2(R,Z,[nbins,nbins],'FaceColor','flat','Normalization','probability');
%     colorbar
%     xAxis = n.XBinEdges;
%     dx = mean( diff(xAxis) );
%     yAxis = n.YBinEdges;
%     dy = mean( diff(yAxis) );
%     axis([min(xAxis) max(xAxis) min(yAxis) max(yAxis)])
%     axis equal
%     view([0 90])
%     colormap(jet(1024))
%     title(['Species: ' num2str(ss) 'Time: ' num2str(ST.time(it))],'Interpreter','latex','FontSize',11)
%     xlabel('$R$','Interpreter','latex','FontSize',16)
%     ylabel('$Z$','Interpreter','latex','FontSize',16)
%     
%     x = ST.params.fields.Ro + ...
%         ST.params.species.ro(ss)*cos(linspace(0,2*pi,100));
%     y = ST.params.species.ro(ss)*sin(linspace(0,2*pi,100));
%     hold on; plot3(x,y,1*ones(size(x)),'k');hold off
%     
%     F(ii) = getframe(gcf);
% end
end

end

function angularMomentum(ST)
Bo = ST.params.fields.Bo;
Ro = ST.params.fields.Ro; % Major radius in meters.
lambda = ST.params.fields.lambda;
qo = ST.params.fields.qo;

st1 = zeros(ST.num_snapshots,ST.params.simulation.num_species);
st2 = zeros(ST.num_snapshots,ST.params.simulation.num_species);
st3 = zeros(ST.num_snapshots,ST.params.simulation.num_species);
st4 = zeros(ST.num_snapshots,ST.params.simulation.num_species);

h = figure('Visible',ST.visible);
set(h,'name','Angular momentum conservation','numbertitle','off')
h1 = figure('Visible',ST.visible);
set(h1,'name','Angular momentum statistics','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    num_part = numel(find(pin==1));
%     num_part = ST.params.species.ppp(ss)*ST.params.simulation.nmpi;
    
    m = ST.params.species.m(ss);
    q = ST.params.species.q(ss);
    
    I = zeros(num_part,ST.num_snapshots);
    err = zeros(1,ST.num_snapshots);
    minerr = zeros(1,ST.num_snapshots);
    maxerr = zeros(1,ST.num_snapshots);
    
    for ii=1:ST.num_snapshots
        X = squeeze( ST.data.(['sp' num2str(ss)]).X(:,pin,ii) );
        V = squeeze( ST.data.(['sp' num2str(ss)]).V(:,pin,ii) );
        gammap = squeeze( ST.data.(['sp' num2str(ss)]).g(pin,ii) )';
        
        
        % Toroidal coordinates
        % r = radius, theta = poloidal angle, phi = toroidal angle
        r = sqrt( (sqrt(X(1,:).^2 + X(2,:).^2) - Ro).^2 + X(3,:).^2 );
        theta = atan2(X(3,:),sqrt(X(1,:).^2 + X(2,:).^2) - Ro);
        theta(theta<0) = theta(theta<0) + 2*pi;
        zeta = atan2(X(1,:),X(2,:));
        zeta(zeta<0) = zeta(zeta<0) + 2*pi;
        % Toroidal coordinates

        eta = r/Ro;
        wc = q*Bo./(m*gammap);

        dzeta_dt = (X(2,:).*V(1,:) - X(1,:).*V(2,:))./( X(1,:).^2 + X(2,:).^2 );
    
        T1 = dzeta_dt.*(1 + eta.*cos(theta)).^2;
        T2 = lambda^2*wc.*log(1 + (r/lambda).^2)/(2*qo*Ro^2);
        I(:,ii) = T1 + T2;
        
  %      disp(I(:,ii))
        
        tmp_vec = (I(:,ii) - I(:,1))./I(:,1);
        err(ii) = mean( tmp_vec );
%         minerr(ii) = min( tmp_vec );
%         maxerr(ii) = max( tmp_vec );
        minerr(ii) = err(ii) - std( tmp_vec );
        maxerr(ii) = err(ii) + std( tmp_vec );
        
        if (~isempty(tmp_vec))
            st1(ii,ss) = mean(tmp_vec);
            st2(ii,ss) = std(tmp_vec);
            %st3(ii,ss) = skewness(tmp_vec);
            %st4(ii,ss) = kurtosis(tmp_vec);
        end
    end
    
    for pp=1:size(I,1)
        I(pp,:) = (I(pp,:) - I(pp,1))/I(pp,1);
    end
    
    figure(h)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    try
        plot(ST.time,err,'k-',ST.time,minerr,'r:',ST.time,maxerr,'r:')
        hold on
        plot(ST.time,I)
        hold off
    catch
    end
    box on
    grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',12)
    ylabel('$\Delta p_\zeta/p_{\zeta 0}$ (\%)','Interpreter','latex','FontSize',12)
    saveas(h,[ST.path 'ang_mom_conservation'],'fig')
    
    figure(h1)
    subplot(4,1,1)
    hold on; plot(ST.time,st1(:,ss)); hold off
    box on; grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$\mu$','Interpreter','latex','FontSize',16)
    subplot(4,1,2)
    hold on; plot(ST.time,st2(:,ss)); hold off
    box on; grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$\sigma$','Interpreter','latex','FontSize',16)
    subplot(4,1,3)
    hold on; plot(ST.time,st3(:,ss)); hold off
    box on; grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$s$','Interpreter','latex','FontSize',16)
    subplot(4,1,4)
    hold on; plot(ST.time,st4(:,ss)); hold off
    box on; grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$k$','Interpreter','latex','FontSize',16)
    saveas(h1,[ST.path 'ang_mom_stat'],'fig')
end


end

function CMF = changeOfMagneticField(ST)
CMF = struct;

CMF.stat.mag = zeros(2,ST.params.simulation.num_species);
CMF.stat.vec = zeros(2,ST.params.simulation.num_species);

for ss=1:ST.params.simulation.num_species   
    q = ST.params.species.q(ss);
    m = ST.params.species.m(ss);
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    aux = find(pin == 1);
    
    Bo = squeeze( sqrt( sum(ST.data.(['sp' num2str(ss)]).B(:,pin,1).^2,1) ) );
    gammap = ST.data.(['sp' num2str(ss)]).g(pin,1)';
    wc = abs(q)*Bo./(gammap*m);
    Tc = 2*pi./wc;
    I = zeros(size(Bo));
    B = zeros(size(Bo));
    DBvec = zeros(size(Bo));
    DBmag = zeros(size(Bo));
    for ii=1:numel(Bo)
        [~,I(ii)] = min( abs(Tc(ii) - ST.time) );
        tmp = squeeze( sqrt( sum(ST.data.(['sp' num2str(ss)]).B(:,aux(ii),1:I(ii)).^2,1) ) );
        DBmag(ii) = max( 100*abs(tmp - Bo(ii))/Bo(ii) );
        
        B1 = squeeze(ST.data.(['sp' num2str(ss)]).B(:,aux(ii),1));
        for jj=1:I(ii)
            B2 = squeeze(ST.data.(['sp' num2str(ss)]).B(:,aux(ii),jj));
            tmp = 100*sqrt(sum((B2-B1).^2))/sqrt(sum(B1.^2));
            if(abs(tmp) > DBvec(ii))
                DBvec(ii) = tmp;
            end
        end
        
%         B(ii) = squeeze( sqrt( sum(ST.data.(['sp' num2str(ss)]).B(:,aux(ii),I(ii)).^2,1) ) );
    end
%     DBmag = 100*abs( B - Bo )./Bo;

    CMF.stat.mag(1,ss) = mean(DBmag);
    CMF.stat.mag(2,ss) = std(DBmag);
    CMF.stat.vec(1,ss) = mean(DBvec);
    CMF.stat.vec(2,ss) = std(DBvec);

    xmag = linspace(min(DBmag),max(DBmag),35);
    xvec = linspace(min(DBvec),max(DBvec),35);
    xAxisMin = min([min(DBmag) min(DBvec)]);
    xAxisMax = max([max(DBmag) max(DBvec)]);

    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,pin,1));
    R = sqrt( sum(X(1:2,:).^2,1) );
    Z = X(3,:);

    S = 12*ones(size(gammap));
    h=figure('Visible',ST.visible,'units','normalized','OuterPosition',[0.1,0.25,0.75,0.4]);
    set(h,'name',['Change of B-field: sp' num2str(ss)],'numbertitle','off')
    subplot(1,3,1)
    histogram(DBvec,xvec)
    hold on
    histogram(DBmag,xmag)
    hold off
    xlim([xAxisMin xAxisMax])
    currentAxis = gca;
    set(currentAxis.Children(:),'FaceAlpha',0.2);
    ylabel('$f(\Delta B/B_0)$','Interpreter','latex','FontSize',16)
    xlabel('$\Delta B/B_0$ ($\%$)','Interpreter','latex','FontSize',16)
    subplot(1,3,2)
    scatter3(R,Z,DBvec,S,DBvec,'square','filled')
    title('$\Delta \mathbf{B}$','Interpreter','latex','FontSize',16)
    colormap(jet(256))
    colorbar
    view([0,90])
    box on; axis on; axis square
    xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    subplot(1,3,3)
    scatter3(R,Z,DBmag,S,DBmag,'square','filled')
    title('$\Delta B$','Interpreter','latex','FontSize',16)
    colormap(jet(256))
    colorbar
    view([0,90])
    box on; axis on; axis square
    xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    saveas(h,[ST.path 'B-field_sp' num2str(ss)],'fig')
%     close(h)
end



end

function RT = radialTransport(ST)
RT = struct;

% cad = ST.params.simulation.output_cadence;
% time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);
Ro = ST.params.fields.Ro;
rc = zeros(1,ST.params.simulation.num_species);

C = colormap(jet(512));
offset = floor(512/ST.params.simulation.num_species);
colour = C(1:offset:end,:);

h1=figure;
set(h1,'name','IC','numbertitle','off')
legends = cell(1,ST.params.simulation.num_species);
for ss=1:ST.params.simulation.num_species
%     h1=figure;
%     set(h1,'name',['Species ' num2str(ss)],'numbertitle','off')
    pin = zeros(1,ST.params.species.ppp(ss)*ST.params.simulation.nmpi);
    for pp=1:ST.params.species.ppp(ss)*ST.params.simulation.nmpi
        X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,pp,:));
        
        % Toroidal coordinates
        % r = radius, theta = poloidal angle, zeta = toroidal angle
        r = sqrt( (sqrt(sum(X(1:2,:).^2,1)) - Ro).^2 + X(3,:).^2 );
        theta = atan2(X(3,:),sqrt(sum(X(1:2,:).^2,1)) - Ro);
        theta(theta<0) = theta(theta<0) + 2*pi;
        theta = 180*theta/pi;
        zeta = atan2(X(1),X(2));
        zeta(zeta<0) = zeta(zeta<0) + 2*pi;
        % Toroidal coordinates
        
        bool = all(ST.data.(['sp' num2str(ss)]).eta(pp,:) < 90);
        
        if all(r < ST.params.fields.a) && bool && all
%             rc(ss) = r(1);
            pin(pp) = 1;
        end
%        figure(h1)
%        hold on
%        plot(theta,r,'.')
%        hold off
    end
    
    t = linspace(0,2*pi,200);
    Rs = ST.params.fields.Ro + ST.params.fields.a*cos(t);
    Zs = ST.params.fields.a*sin(t);

    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,logical(pin),1));
    R = sqrt( sum(X(1:2,:).^2,1) );
    Z = X(3,:);
    
    figure(h1)
%     subplot(2,1,1)
    hold on
    plot(R,Z,'k.','MarkerSize',12,'MarkerFaceColor',colour(ss,:),'MarkerEdgeColor',colour(ss,:))
    hold off
%     hold on
%     plot(Rs,Zs,'r')
%     box on
%     axis on
%     axis equal
%     xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
%     ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)

    
%     X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,logical(pin),end));
%     R = sqrt( sum(X(1:2,:).^2,1) );
%     Z = X(3,:);
%     
%     subplot(2,1,2)
%     plot(R,Z,'b.',Rs,Zs,'r')
%     box on
%     axis on
%     axis equal
%     xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
%     ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    legends{ss} = ['$\eta_0 =$' num2str(ST.params.species.etao(ss)) '$^\circ$'];
end

figure(h1)
% subplot(2,1,1)
legend(legends)
hold on
plot(Rs,Zs,'r')
box on
axis on
axis square
xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)

RT.rc = rc;

end

function CP = confined_particles(ST)
CP = struct;

tmax = max(ST.time);
tmin = min(ST.time);

if isfield(ST.params.fields,'a')
    t = linspace(0,2*pi,200);
    Rs = ST.params.fields.Ro + ST.params.fields.a*cos(t);
    Zs = ST.params.fields.a*sin(t);
end

CP.confined = zeros(1,ST.params.simulation.num_species);

h0 = figure('Visible',ST.visible);
set(h0,'name','Particle loss','numbertitle','off');

C = colormap(h0,jet(512));
offset = floor(512/ST.params.simulation.num_species);
colour = C(1:offset:end,:);

legends = cell(1,ST.params.simulation.num_species);
for ss=1:ST.params.simulation.num_species
    numConfPart = sum(ST.data.(['sp' num2str(ss)]).flag(:,1),1);
    confinedParticles = ...
        100*sum(ST.data.(['sp' num2str(ss)]).flag,1)/numConfPart;
    figure(h0)
    hold on
    plot(ST.time,confinedParticles)
    hold off
    
    legends{ss} = ['$\eta_0 =$' num2str(ST.params.species.etao(ss)) '$^\circ$'];
    
    CP.confined(ss) = confinedParticles(end);
end


figure(h0)
legend(legends,'interpreter','latex','FontSize',12)
box on
axis on
axis square
xlabel('Time $t$ (sec)','Interpreter','latex','FontSize',16)
ylabel('$\%$ of confined RE','Interpreter','latex','FontSize',16)
saveas(h0,[ST.path 'particle_loss'],'fig')
% close(h0)

h1=figure;
set(h1,'Visible',ST.visible,'name','IC','numbertitle','off')
legends = cell(1,ST.params.simulation.num_species);
for ss=1:ST.params.simulation.num_species   
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    bool = pin & passing;
%     bool = pin;
    
    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,bool,1));
    R = sqrt( sum(X(1:2,:).^2,1) );
    Z = X(3,:);

    figure(h1)
    subplot(1,2,1)
    hold on
    plot(R,Z,'.','MarkerSize',8,'MarkerFaceColor',colour(ss,:),'MarkerEdgeColor',colour(ss,:))
%     his=histogram2(R,Z,'FaceColor','flat',...
%         'DisplayStyle','tile','ShowEmptyBins','on','LineStyle','none');
%     colormap(jet(1024));caxis([min(min(his.Values)) max(max(his.Values))])
%     colorbar
    hold off
    legends{ss} = ['$\eta_0 =$' num2str(ST.params.species.etao(ss)) '$^\circ$'];
end

figure(h1)
subplot(1,2,1)
legend(legends,'interpreter','latex','FontSize',12)
if isfield(ST.params.fields,'a')
    hold on
    plot(Rs,Zs,'r')
    hold off
end
box on
axis on
axis square
xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)

legends = cell(1,ST.params.simulation.num_species);
for ss=ST.params.simulation.num_species:-1:1
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    trapped = logical( any(ST.data.(['sp' num2str(ss)]).eta > 90,2) );
    bool = pin & trapped;
    
    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,bool,1));
    R = sqrt( sum(X(1:2,:).^2,1) );
    Z = X(3,:);
    
    figure(h1)
    subplot(1,2,2)
    hold on
    plot(R,Z,'.','MarkerSize',8,'MarkerFaceColor',colour(ss,:),'MarkerEdgeColor',colour(ss,:))
%     plot(R,Z,'.','MarkerSize',10,'MarkerFaceColor',colour(ss,:),'MarkerEdgeColor',colour(ss,:))
%     his=histogram2(R,Z,'FaceColor','flat',...
%         'DisplayStyle','tile','ShowEmptyBins','on','LineStyle','none');
%     colormap(jet(1024));caxis([min(min(his.Values)) max(max(his.Values))])
%     colorbar
    hold off
    legends{ST.params.simulation.num_species + 1 -ss} = ['$\eta_0 =$' num2str(ST.params.species.etao(ss)) '$^\circ$'];
end


figure(h1)
subplot(1,2,2)
legend(legends,'interpreter','latex','FontSize',12)
if isfield(ST.params.fields,'a')
    hold on
    plot(Rs,Zs,'r')
    box on
end
axis on
axis square
xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
saveas(h1,[ST.path 'IC'],'fig')
% close(h1)
end

function energyLimit(ST)
% cad = ST.params.simulation.output_cadence;
% time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);
tmax = max(ST.time);
tmin = min(ST.time);

h1=figure;
set(h1,'name','Energy vs. pitch angle','numbertitle','off')
h2=figure;
set(h2,'name','Energy vs. time','numbertitle','off')
legends = cell(1,ST.params.simulation.num_species);

for ss=1:ST.params.simulation.num_species   
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    
    energy = mean(ST.data.(['sp' num2str(ss)]).g(pin,:),1);
    energy = ST.params.species.m(ss)*ST.params.scales.v^2*energy/ST.params.scales.q;
    energy = energy/1E6;
    pitch = mean(ST.data.(['sp' num2str(ss)]).eta(pin,:),1);
    
    figure(h1)
    hold on
%     plot(pitch(1),energy(1),'ko',pitch,energy)
    plot(pitch,energy)
    hold off
        
    figure(h2)
    hold on
    plot(ST.time,energy)
    hold off
    legends{ss} = ['$\eta_0 =$' num2str(ST.params.species.etao(ss)) '$^\circ$'];
end

figure(h1)
legend(legends,'interpreter','latex','FontSize',12)
grid on; box on; axis on; axis square
xlabel('Pitch angle $\eta$ ($^\circ$)','Interpreter','latex','FontSize',16)
ylabel('$\mathcal{E}_0$ (MeV)','Interpreter','latex','FontSize',16)
saveas(h1,[ST.path 'Energy_vs_pitch'],'fig')


figure(h2)
legend(legends,'interpreter','latex','FontSize',12)
grid on; box on; axis on; axis square
xlabel('Time $t$ (sec)','Interpreter','latex','FontSize',16)
ylabel('$\mathcal{E}_0$ (MeV)','Interpreter','latex','FontSize',16)
saveas(h2,[ST.path 'mean_energy_vs_time'],'fig')
end

function B = analyticalB(ST,X)
% Analytical magnetic field
% X is a vector X(1)=x, X(2)=y, X(3)=z.

narginchk(1,2);

% Parameters of the analytical magnetic field
Bo = ST.params.fields.Bo;
a = ST.params.fields.a;
Ro = ST.params.fields.Ro;
qa = ST.params.fields.qa;
lamb = ST.params.fields.lambda;
try
    co = ST.params.fields.co;
catch
    qo = ST.params.fields.qa;
    co = a/lamb;
end
Bpo = ST.params.fields.Bpo;
% Parameters of the analytical magnetic field


% Toroidal coordinates
% r = radius, theta = poloidal angle, zeta = toroidal angle
r = squeeze(sqrt( (sqrt(sum(X(1:2,:,:).^2,1)) - Ro).^2 + X(3,:,:).^2));
theta = atan2(squeeze(X(3,:,:)),squeeze(sqrt(sum(X(1:2,:,:).^2,1)) - Ro));
theta(theta<0) = theta(theta<0) + 2*pi;
zeta = atan2(squeeze(X(1,:,:)),squeeze(X(2,:,:)));
zeta(zeta<0) = zeta(zeta<0) + 2*pi;
% Toroidal coordinates

% Poloidal magnetic field
% Minus sign = TEXTOR
% Plus sign = default
Bp = Bpo*(r/lamb)./( 1 + (r/lamb).^2 );

eta = r/Ro;
Br = 1./( 1 + eta.*cos(theta) );

Bx = Br.*( Bo*cos(zeta) - Bp.*sin(theta).*sin(zeta) );
By = -Br.*( Bo*sin(zeta) + Bp.*sin(theta).*cos(zeta) );
Bz = Br.*Bp.*cos(theta);

B = zeros(3,size(r,1),size(r,2));
B(1,:,:) = Bx;
B(2,:,:) = By;
B(3,:,:) = Bz;

end

function E = analyticalE(ST,X)
% Analytical magnetic field
% X is a vector X(1)=x, X(2)=y, X(3)=z.

narginchk(1,2);

% Parameters of the analytical magnetic field
Eo = ST.params.fields.Eo;
Ro = ST.params.fields.Ro;
% Parameters of the analytical magnetic field


% Toroidal coordinates
% r = radius, theta = poloidal angle, zeta = toroidal angle
r = squeeze(sqrt( (sqrt(sum(X(1:2,:,:).^2,1)) - Ro).^2 + X(3,:,:).^2));
theta = atan2(squeeze(X(3,:,:)),squeeze(sqrt(sum(X(1:2,:,:).^2,1)) - Ro));
theta(theta<0) = theta(theta<0) + 2*pi;
zeta = atan2(squeeze(X(1,:,:)),squeeze(X(2,:,:)));
zeta(zeta<0) = zeta(zeta<0) + 2*pi;
% Toroidal coordinates

% Poloidal magnetic field
eta = r/Ro;
Ezeta = Eo./( 1 + eta.*cos(theta) );

Ex = Ezeta.*cos(zeta);
Ey = -Ezeta.*sin(zeta);
Ez = 0;

E = zeros(3,size(r,1),size(r,2));
E(1,:,:) = Ex;
E(2,:,:) = Ey;
E(3,:,:) = Ez;

end

function PR = LarmorVsLL(ST)
PR = struct;

kB = 1.38E-23; % Boltzmann constant
Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
mu0 = (4E-7)*pi; % Magnetic permeability
ep = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass
% % % % % % % % % % %

PR.models = zeros(3,ST.params.simulation.num_species);

h = figure('Visible',ST.visible);
set(h,'name','Model comparison: ratios','numbertitle','off')
h1 = figure('Visible',ST.visible);
set(h1,'name','Scatter plot: comparison','numbertitle','off')
h2 = figure('Visible',ST.visible);
set(h2,'name','Scatter plot: actual radiation','numbertitle','off')
h3 = figure('Visible',ST.visible);
set(h3,'name','Radiation bar plots','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    q = abs(ST.params.species.q(ss));
    m = ST.params.species.m(ss);
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    bool = pin & passing;
    aux = find(bool == 1);
    S = numel(aux);
    
    V = ST.data.(['sp' num2str(ss)]).V(:,bool,:);
    v = squeeze( sqrt( sum(V.^2,1) ) );
    gammap = ST.data.(['sp' num2str(ss)]).g(bool,:);
    eta = pi*ST.data.(['sp' num2str(ss)]).eta(bool,:)/180;
    
    gammapo = repmat(ST.params.species.go(ss),S,ST.num_snapshots);
    vo = ST.params.scales.v*sqrt(1 - 1./gammapo.^2);
    etao = pi*repmat(ST.params.species.etao(ss),S,ST.num_snapshots)/180;
    
%     gammapo = gammap;
%     etao = eta;
%     vo = v;
    
    X = ST.data.(['sp' num2str(ss)]).X(:,bool,:);
    try
        B = ST.data.(['sp' num2str(ss)]).B(:,bool,:);
    catch
        B = analyticalB(ST,X);
    end
    E = analyticalE(ST,X);
    
    vec_mag = zeros(size(gammap));
    for it=1:size(X,3)
        for ii=1:size(gammap,1)
            VxE = cross(squeeze(V(:,ii,it)),squeeze(E(:,ii,it)));
            VxB = cross(squeeze(V(:,ii,it)),squeeze(B(:,ii,it)));
            VxVxB = cross(squeeze(V(:,ii,it)),VxB);
            vec = VxE + VxVxB;
            vec_mag(ii,it) = sqrt( vec'*vec );
        end
    end

    % Approximation of <1/R^2> ~ sin^4(eta)/rg^2
    vperp = vo.*sin(etao);
    wc = q*ST.params.fields.Bo./(gammapo*m);
    rg = vperp./wc;
    kappa2 = sin(etao).^4./rg.^2;
 
    % Actual curvature
    kappa = q*vec_mag./(m*gammap.*v.^3);    

    % Landau-Lifshiftz radiation formula
    PR_LL = abs(ST.data.(['sp' num2str(ss)]).Prad(bool,:));
    
    % Larmor approximation using actual curvature
%     PR_L = 2*Kc*q^2*(gammap.*v).^4.*kappa.^2/(3*c^3);
    Tr = 6*pi*ep*(m*ST.params.scales.v)^3./(q^4*squeeze(sum(B.^2,1)));
    PR_L = gammap.*v.*(m*gammap.*v).*sin(eta).^2./Tr;
    
    % Larmor approximation using approximation for curvature
%     PR_app = 2*Kc*q^2*(gammapo.*v).^4.*kappa2/(3*c^3);
    Tr = 6*pi*ep*(m*ST.params.scales.v)^3/(q^4*ST.params.fields.Bo^2);
    PR_app = gammapo.*vo.*(m*gammapo.*vo).*sin(etao).^2/Tr;
    
    figure(h1)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    hold on
    plot(squeeze(180*eta(:,end)/pi),squeeze(PR_L(:,end)),'b.','MarkerSize',2)
    plot(squeeze(180*etao(:,end)/pi),squeeze(PR_app(:,end)),'r.','MarkerSize',2)
%     plot(squeeze(180*eta(:,end)/pi),squeeze(PR_app(:,end)),'r.','MarkerSize',2)
    hold off
    box on; grid on;
    xlabel('$\theta$ ($^\circ$)','Interpreter','latex','FontSize',16)
    ylabel('$P_R$','Interpreter','latex','FontSize',16)
    
    figure(h2)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    plot(squeeze(180*eta(:,end)/pi),squeeze(PR_LL(:,end)),'k.','MarkerSize',4)
    box on; grid on;
    xlabel('$\theta$ ($^\circ$)','Interpreter','latex','FontSize',16)
    ylabel('$P_R$','Interpreter','latex','FontSize',16)
    
    try
        minVal = min( eta(:,end) );
        maxVal = max( eta(:,end) );
        x = linspace(minVal,maxVal,20);
        Dx = mean(diff(x));
        
        [fx,~] = hist(eta(:,end),x);
%         fx = fx/(Dx*sum(fx));
    catch
    end
    
    
    try
        minVal = min( PR_LL(:,end) );
        maxVal = max( PR_LL(:,end) );
        z = linspace(minVal,maxVal,20);
        Dz = mean(diff(z));
        
        [fz,~] = hist(PR_LL(:,end),z);
%         fz = fz/(Dz*sum(fz));
    catch
    end
    
    figure(h3)
    offset = 2*(double(ss) - 1);
    subplot(double(ST.params.simulation.num_species),2,offset + 1)
    bar(180*x/pi,fx)
    xlabel('$\theta$ ($^\circ$)','Interpreter','latex','FontSize',16)
    subplot(double(ST.params.simulation.num_species),2,offset + 2)
    bar(z,fz)
    xlabel('$P_R$ (Watts/electron)','Interpreter','latex','FontSize',16)
    
%     PR_LL = mean(PR_LL,1);
    PR_LL = sum(PR_LL,1);
    PR_L = mean(PR_L,1);
    PR_app = mean(PR_app,1);
    
    RATIO1 = PR_L./PR_LL;
    RATIO2 = PR_app./PR_LL;   
    PR.models(:,ss) = [PR_LL(end), PR_L(end), PR_app(end)];
    
    figure(h)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    plot(ST.time,RATIO1,'k',ST.time,RATIO2,'r')
    box on; grid on;
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$P_R$','Interpreter','latex','FontSize',12)
end
saveas(h,[ST.path 'radiation_ratios'],'fig')
saveas(h1,[ST.path 'scatter_plot_comparison'],'fig')
saveas(h2,[ST.path 'scatter_plot_actual_radiation'],'fig')
saveas(h3,[ST.path 'radiation_pdfs'],'fig')
end

function stackedPlots(ST,nbins)
tmax = max(ST.time);
tmin = min(ST.time);

fx = zeros(ST.params.simulation.num_species,nbins,ST.num_snapshots);
Prad = zeros(ST.params.simulation.num_species,ST.num_snapshots);
x = zeros(ST.params.simulation.num_species,nbins);

h2 = figure;
set(h2,'name','PDF of radiated power','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
%     data = ST.data.(['sp' num2str(ss)]).eta(pin,:);
    data = abs(ST.data.(['sp' num2str(ss)]).Prad(pin,:));

    minVal = min(min( data ));
    maxVal = max(max( data ));
    x(ss,:) = linspace(minVal,maxVal,nbins);
    
    for ii=1:ST.num_snapshots
        [fx(ss,:,ii),~] = hist(data(:,ii),x(ss,:));
        Prad(ss,ii) = mean(diff(x(ss,:)))*sum(x(ss,:).*fx(ss,:,ii));
    end
    
    figure(h2)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
%     surf(ST.time,squeeze(x(ss,:)),log10(squeeze(fx(ss,:,:))),'LineStyle','none')
    plot(ST.time,Prad(ss,:))
    axis([tmin tmax min(Prad(ss,:)) max(Prad(ss,:))])
%     axis([tmin tmax minVal maxVal])
    box on
    axis on
    xlabel('Time $t$ (sec)','Interpreter','latex','FontSize',16)
    ylabel('$P_{rad}$ (Watts/electron)','Interpreter','latex','FontSize',16)
    colormap(jet(1024))
    colorbar
end


end

function scatterPlots(ST)
kB = 1.38E-23; % Boltzmann constant
Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
mu0 = (4E-7)*pi; % Magnetic permeability
ep = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass
% % % % % % % % % % %

h=figure;
set(h,'name','Scatter: E vs pitch','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    q = abs(ST.params.species.q(ss));
    m = ST.params.species.m(ss);
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    eta = ST.data.(['sp' num2str(ss)]).eta(pin,:);
    Et = ST.data.(['sp' num2str(ss)]).g(pin,:)*m*ST.params.scales.v^2/(abs(q)*1E6);
    Eo = ST.params.species.Eo(ss)/1E6;
    etao = ST.params.species.etao(ss);

    eta_mean = mean(eta,1);
    Et_mean = mean(Et,1);

    figure(h)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    plot(squeeze(eta(:,end)),squeeze(Et(:,end)),'.','MarkerFaceColor',[0,0.45,0.74],...
        'MarkerEdgeColor',[0,0.45,0.74])
%     hold on
%     plot(eta_mean,Et_mean,'k',...
%         etao,Eo,'rs','MarkerFaceColor','r','MarkerSize',10)
%     hold off
    box on; grid on;
    xlabel('$\eta$ ($^\circ$)','Interpreter','latex','FontSize',16)
    ylabel('$\mathcal{E}$ (MeV)','Interpreter','latex','FontSize',16)
end
saveas(h,[ST.path 'scatter_plot_E_vs_pitch'],'fig')

end

function P = synchrotronSpectrum(ST,filtered)
disp('Calculating spectrum of synchrotron radiation...')
P = struct;

geometry = 'cylindrical';

upper_integration_limit = 200.0;

N = 250;
lambda_min = 1E-9;% in meters
lambda_max = 5000E-9;% in meters
% 
% lambda_min = 450E-9;% in meters
% lambda_max = 950E-9;% in meters
% lambda_min = 907E-9;% in meters
% lambda_max = 917E-9;% in meters
% lambda_min = 742E-9;% in meters
% lambda_max = 752E-9;% in meters

% lambda_min = 10E-9;% in meters
% lambda_max = 2E-5;% in meters

lambda_camera = linspace(lambda_min,lambda_max,N);
Dlambda_camera = mean(diff(lambda_camera));

rmin = 0;
try
    rmax = ST.params.fields.a;
    Nr = 25;
    Ntheta = 80;
    
    Rmin = 0.9;
    Rmax = 2.1;
    
    Zmin = -0.6;
    Zmax = 0.6;
    
    NR = 50;
    NZ = 50;
catch
    %     Ro = ST.params.fields.Ro;
    %     rmax = max([max(ST.params.fields.R) - Ro, Ro - min(ST.params.fields.R)]);
    rmax = 1.2;
    Nr = 40;
    Ntheta = 80;
    
    Rmin = 0.9;
    Rmax = 2.6;
    
    Zmin = -2.0;
    Zmax = 2.0;
    
    NR = 100;%100
    NZ = 235;%235
end

if strcmp(geometry,'poloidal')
    Psyn_total = zeros(Ntheta,Nr);
else
    Psyn_total = zeros(NZ,NR);
end


lambda_camera = 1E2*lambda_camera;

num_species = double(ST.params.simulation.num_species);

fh = figure;
spectrum_figure = figure;
numPanels = ceil(sqrt(num_species + 1));

L = cell(1,num_species);

% Poloidal distribution of the total radiated power
for ss=1:num_species
% for ss=6:6
    q = abs(ST.params.species.q(ss));
    m = ST.params.species.m(ss);
    Ro = ST.params.fields.Ro;
    
    if strcmp(geometry,'poloidal')
        Psyn = zeros(Ntheta,Nr);
    else
        Psyn = zeros(NZ,NR);
    end
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    bool = pin;% & passing;
    
    X = [];
    V = [];
    gammap = [];
    Prad = [];
    eta = [];
    for ii=1:ST.num_snapshots
        X = [X,ST.data.(['sp' num2str(ss)]).X(:,bool,ii)];
        V = [V,ST.data.(['sp' num2str(ss)]).V(:,bool,ii)];
        gammap = [gammap;ST.data.(['sp' num2str(ss)]).g(bool,ii)];
        eta = [eta;pi*ST.data.(['sp' num2str(ss)]).eta(bool,ii)/180];
        Prad = [Prad;abs(ST.data.(['sp' num2str(ss)]).Prad(bool,ii))];
    end
    v = squeeze( sqrt( sum(V.^2,1) ) )';
    
    [vp,~,camPos] = identifyVisibleParticles(X,V,gammap,false);
%     vp = true(size(gammap));
%     psi = ones(size(gammap));
    
    X(:,~vp) = [];
    V(:,~vp) = [];
    gammap(~vp) = [];
    v(~vp) = [];
    eta(~vp) = [];
    Prad(~vp) = [];
    
    numPart = size(camPos,1);
    
    if strcmp(geometry,'poloidal')
        % Toroidal coordinates
        % r = radius, theta = poloidal angle, zeta = toroidal angle
        r = squeeze(sqrt( (sqrt(sum(X(1:2,:).^2,1)) - Ro).^2 + X(3,:).^2));
        theta = atan2(squeeze(X(3,:)),squeeze(sqrt(sum(X(1:2,:).^2,1)) - Ro));
        theta(theta<0) = theta(theta<0) + 2*pi;

        Dr = (rmax - rmin)/Nr;
        r_grid = 0.5*Dr + (0:1:(Nr-1))*Dr;

        Dtheta = 2*pi/Ntheta;
        theta_grid = 0.5*Dtheta + (0:1:(Ntheta-1))*Dtheta;

        ir = floor(r/Dr) + 1;
        itheta = floor(theta/Dtheta) + 1;

        % % % Set-up of the grid of the poloidal plane
        x_grid = zeros(Ntheta,Nr);
        y_grid = zeros(Ntheta,Nr);
        for ii=1:Nr
            for jj=1:Ntheta
                x_grid(jj,ii) = Ro + r_grid(ii)*cos(theta_grid(jj));
                y_grid(jj,ii) = r_grid(ii)*sin(theta_grid(jj));
            end
        end
        % Toroidal coordinates
    else
        % Cylindrical coordinates
        R = sqrt(sum(X(1:2,:).^2,1));
        Z = X(3,:);

        DR = (Rmax - Rmin)/NR;
        DZ = (Zmax - Zmin)/NZ;

        R_grid = Rmin + 0.5*DR + (0:1:(NR-1))*DR;
        Z_grid = Zmin + 0.5*DZ + (0:1:(NZ-1))*DZ;

        iR = floor((R - Rmin)/DR) + 1;
        iZ = floor((Z - Zmin)/DZ) + 1;
        % Cylindrical coordinates
    end
    
    % % % Option for calculating the total Psyn WITHOUT wavelength filtering
    if (~filtered)
        for ii=1:numPart
            try
                if strcmp(geometry,'poloidal')
                    Psyn(itheta(ii),ir(ii)) = Psyn(itheta(ii),ir(ii)) + ...
                        Prad(ii);
                    
                    Psyn_total(itheta(ii),ir(ii)) = Psyn_total(itheta(ii),ir(ii)) + ...
                        Prad(ii);
                else
                    Psyn(iZ(ii),iR(ii)) = ...
                        Psyn(iZ(ii),iR(ii)) + Prad(ii);
                    
                    Psyn_total(iZ(ii),iR(ii)) = ...
                        Psyn_total(iZ(ii),iR(ii)) + Prad(ii);
                end
            catch
                disp('An issue calculating poloidal plane')
            end
        end
        
        figure(fh)
        subplot(numPanels,numPanels,ss)
        if strcmp(geometry,'poloidal')
            surf(x_grid,y_grid,Psyn,'LineStyle','none')
            axis([min(x_grid) max(x_grid) min(y_grid) max(y_grid)])
        else
            contourf(R_grid,Z_grid,Psyn,'LineStyle','none')
            axis([min(R_grid) max(R_grid) min(Z_grid) max(Z_grid)])
        end
        colormap(jet(512))
        h = colorbar;
        ylabel(h,'$P_{syn}$ (Watts)','Interpreter','latex','FontSize',16)
        view([0,90])
        axis equal; box on
        shading interp
        xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
        ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
        title(['$\theta_0 = $' num2str(ST.params.species.etao(ss)) '$^\circ$'],...
            'Interpreter','latex','FontSize',16)
        
    else % WITH filter
        
        try
            B = ST.data.(['sp' num2str(ss)]).B(:,bool,it);
            E = zeros(size(B));
        catch
            B = analyticalB(ST,X);
            E = analyticalE(ST,X);
        end

        vec_mag = zeros(numPart,1);
        Binormal = zeros(numPart,3);
        
        for ii=1:numPart
            VxE = cross(squeeze(V(:,ii)),squeeze(E(:,ii)));
            VxB = cross(squeeze(V(:,ii)),squeeze(B(:,ii)));
            VxVxB = cross(squeeze(V(:,ii)),VxB);
            vec = VxE + VxVxB;
            vec_mag(ii) = sqrt( vec'*vec );
            
            Binormal(ii,:) = vec/vec_mag(ii);
        end
        
        psi = abs(pi/2 - acos(dot(camPos',Binormal')));
        
        % Actual curvature
        k = q*vec_mag./(m*gammap.*v.^3);

        % % % % NO CUTOFF % % %
%         k = ...
%             q*ST.params.fields.Bo*sin(pi*ST.params.species.etao(ss)/180)...
%             /(ST.params.species.go(ss)*ST.params.species.m(ss)*v(1))*ones(size(v));
        
        % % % % Beyond this point all variables are in cgs units % % % %
        c = 1E2*ST.params.scales.v;
        qe = 3E9*q;
        m = 1E3*m;
%         k = k/1E2;
        
        lambdac = (4/3)*pi*(1./gammap).^3./k;
        
%         lambda_camera = linspace(1E-7,5*lambdac(1),N); % % % % NO CUTOFF % % %
        
        I = find(lambdac > lambda_min);
        numEmittingPart = numel(I);
         
        Psyn_camera = zeros(numEmittingPart,N);
        disp('Decomposing radiation in wavelengths...')
        
        C0 = 4*pi*c*qe^2/sqrt(3);
        fun = @(x) besselk(5/3,x);
        y = (gammap.*psi').^2;
        for ii=1:numEmittingPart
            ind = I(ii);
            for jj=1:N
                lower_integration_limit = lambdac(ind)/lambda_camera(jj);
%                 if (lambda_camera(jj) < lambdac(ind)) && (lower_integration_limit < upper_integration_limit)             
%                 if (lower_integration_limit < upper_integration_limit)
%                     Q = integral(fun,lower_integration_limit,upper_integration_limit);
%                     C1 = 1/(gammap(ind)^2*lambda_camera(jj)^3);
%                     Psyn_camera(ii,jj) =  C0*C1*Q;
%                 end

                if (lambda_camera(jj) < lambdac(ind)) && isfinite(lambdac(ind))
                    zeta = 0.5*lower_integration_limit*(1 + y(ind))^(3/2);
                    D0 = 3*c*qe^2*k(ind)/(2*pi*lambda_camera(jj)^2);
                    
                    Psyn_camera(ii,jj) = ...
                        D0*lower_integration_limit^2*gammap(ind)^2*(1 + y(ind))^2*(besselk(2/3,zeta).^2 + ...
                        (y(ind)/(1 + y(ind)))*besselk(1/3,zeta).^2);
                end
            end
        end
        Psyn_camera = Psyn_camera/numEmittingPart;
                
        % With the data of Psyn_camera we make the poloidal plane
        % 'photograph'.
        disp('Calculating poloidal plane...')
        for ii=1:numEmittingPart
            ind = I(ii);
            try
                Psyn_integrated = trapz(lambda_camera,Psyn_camera(ii,:));
                if strcmp(geometry,'poloidal')
                    Psyn(itheta(ind),ir(ind)) = ...
                        Psyn(itheta(ind),ir(ind)) + Psyn_integrated;
                    
                    Psyn_total(itheta(ind),ir(ind)) = ...
                        Psyn_total(itheta(ind),ir(ind)) + Psyn_integrated;
                else
                    Psyn(iZ(ind),iR(ind)) = ...
                        Psyn(iZ(ind),iR(ind)) + Psyn_integrated;
                    
                    Psyn_total(iZ(ind),iR(ind)) = ...
                        Psyn_total(iZ(ind),iR(ind)) + Psyn_integrated;
                end
            catch
                disp('An issue calculating poloidal plane')
            end
        end
        
        lch = 1E7;
        Pch = 1E-7;
        
        lambdac = lch*lambdac;
        Psyn_camera = Pch*Psyn_camera;
        Psyn = Pch*Psyn;
        k = 1E2*k;
        
        % % % % Beyond this point all variables are in SI units % % % %      

        figure(fh)
        subplot(numPanels,numPanels,ss)
        if strcmp(geometry,'poloidal')
            surf(x_grid,y_grid,Psyn,'LineStyle','none')
%             axis([min(x_grid) max(x_grid) min(y_grid) max(y_grid)])
        else
%             surfc(R_grid,Z_grid,Psyn,'LineStyle','none')
            contourf(R_grid,Z_grid,Psyn,'LineStyle','none')
%             axis([min(R_grid) max(R_grid) min(Z_grid) max(Z_grid)])
        end
        colormap(jet(512))
        h = colorbar;
        ylabel(h,'$P_{syn}$ (Watts)','Interpreter','latex','FontSize',16)
        view([0,90])
        axis equal; box on
        shading interp
        xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
        ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
        title(['$\theta_0 = $' num2str(ST.params.species.etao(ss)) '$^\circ$'],...
            'Interpreter','latex','FontSize',16)
        
        figure(spectrum_figure)
        hold on
        plot(lch*lambda_camera,sum(Psyn_camera,1)/lch)
        hold off
        grid on;box on
        xlabel('$\lambda$ (nm)','Interpreter','latex','FontSize',16)
        ylabel('$P_{syn}$ (Watts/m)','Interpreter','latex','FontSize',16)
        L{ss} = ['$\theta_0 = $' num2str(ST.params.species.etao(ss)) '$^\circ$'];
    end
    
end



% % % % Final figures % % % %
if (~filtered)
    figure(fh)
    subplot(numPanels,numPanels,num_species+1)
    if strcmp(geometry,'poloidal')
        surf(x_grid,y_grid,Psyn_total_pol,'LineStyle','none')
%         axis([min(x_grid) max(x_grid) min(y_grid) max(y_grid)])
    else
        surf(R_grid,Z_grid,Psyn_total,'LineStyle','none')
%         axis([min(R_grid) max(R_grid) min(Z_grid) max(Z_grid)])
    end
    colormap(jet(512))
    h = colorbar;
    ylabel(h,'$P_{syn}$ (Watts)','Interpreter','latex','FontSize',16)
    view([0,90])
    axis equal; box on
    shading interp
    xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    title('Total $P_{syn}$','Interpreter','latex','FontSize',16)
else
    Psyn_total = 1E-7*Psyn_total;
    
    figure(fh)
    subplot(numPanels,numPanels,num_species+1)
    if strcmp(geometry,'poloidal')
        surf(x_grid,y_grid,Psyn_total_pol,'LineStyle','none')
        axis([min(x_grid) max(x_grid) min(y_grid) max(y_grid)])
    else
        surf(R_grid,Z_grid,Psyn_total,'LineStyle','none')
        axis([min(R_grid) max(R_grid) min(Z_grid) max(Z_grid)])
    end
    colormap(jet(512))
    h = colorbar;
    ylabel(h,'$P_{syn}$ (Watts)','Interpreter','latex','FontSize',16)
    view([0,90])
    axis equal; box on
    shading interp
    xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    title('Total $P_{syn}$','Interpreter','latex','FontSize',16)
    
    figure(spectrum_figure)
    legend(L,'Interpreter','latex','FontSize',12)
end

disp('Spectrum of synchrotron radiation: done!')
end

function [vp,psi,camPos] = identifyVisibleParticles(X,V,gammap,option)

% Radial position of inner wall
Riw = 1; % in meters

% Radial and vertical position of the camera
Rc = 2.38; % in meters
Zc = 0.0; % in meters
% Zc = 0;

np = numel(gammap);

% mea = maximum emission angle of synchrotron radiation
mea = 1./gammap; % in rad. \psi ~ 1/\gamma

xo = X(1,:);
yo = X(2,:);
zo = X(3,:);

Ro = sqrt(sum(X(1:2,:).^2,1));

v = sqrt(sum(V.^2,1));
vox = V(1,:)./v;
voy = V(2,:)./v;
voz = V(3,:)./v;

tmp_cam = zeros(np,3);

% First we find the Z position where V hits the wall at Rc
theta_f = zeros(1,np);
% Z_f = zeros(1,np);
hitInnerWall = false(1,np);
for ii=1:np
    if (Ro(ii) < Rc)
        p = zeros(1,3);
        
        % polinomial coefficients p(1)*x^2 + p(2)*x + p(3) = 0
        p(1) = vox(ii)^2 + voy(ii)^2;
        p(2) = 2*(xo(ii)*vox(ii) + yo(ii)*voy(ii));
        p(3) = xo(ii)^2 + yo(ii)^2 - Rc^2;
        
        r = roots(p);
        if all(r>0) || all(r<0)
            disp(['Something wrong at ii=' num2str(ii)])
            hitInnerWall(ii) = true;
        else
            t = max(r);
            
            % Cartesian coordinates where the tangent vector intersects the
            % outter wall of the axisymmetric device
            X_f = xo(ii) + vox(ii)*t;
            Y_f = yo(ii) + voy(ii)*t;
%             Z_f = zo(ii) + voz(ii)*t;
            
            theta_f(ii) = atan2(Y_f,X_f);
            if (theta_f(ii) < 0)
                theta_f(ii) = theta_f(ii) + 2*pi;
            end
            
            X_cam = Rc*cos(theta_f(ii));
            Y_cam = Rc*sin(theta_f(ii));
            
            tmp_cam(ii,:) = [X_cam - xo(ii), Y_cam - yo(ii), Zc - zo(ii)];
            tmp_cam(ii,:) = tmp_cam(ii,:)/sqrt(tmp_cam(ii,:)*tmp_cam(ii,:)');
            
            % % % Check if the unitary velocity vector hits the inner wall
            p(3) = xo(ii)^2 + yo(ii)^2 - Riw^2; 
            r = roots(p);
            if isreal(r) && any(r>0)
                hitInnerWall(ii) = true;
            end
        end
    else
        % The particle hits the outer wall
        hitInnerWall(ii) = true;
    end
end

I = find(hitInnerWall == false);

% Then, we calculate the angle between V and the position of the camera
xc = Rc*cos(theta_f(I));
yc = Rc*sin(theta_f(I));

ax = xc - xo(I);
ay = yc - yo(I);
az = Zc - zo(I);

a = sqrt(ax.^2 + ay.^2 + az.^2);

ax = ax./a;
ay = ay./a;
az = az./a;

psi = acos(ax.*vox(I) + ay.*voy(I) + az.*voz(I))';
visible = psi <= mea(I);

if ( option )
    Ro = sqrt(xo(I(visible)).^2 + yo(I(visible)).^2);
    figure
    plot(Ro,zo(I(visible)),'g.','MarkerSize',18)
end

vp = false(1,np);
vp(I(visible)) = true;

psi(~visible) = [];

camPos = tmp_cam(I(visible),:);

end

function [ip,theta_f] = findVisibleParticles(X,V,angle,camera_params,option)
disp('Finding visible particles...');
% Radial and vertical position of the camera
Rc = camera_params.position(1); % in meters
Zc = camera_params.position(2); % in meters

np = numel(angle);

xo = X(1,:);
yo = X(2,:);
zo = X(3,:);

v = sqrt(sum(V.^2,1));
vox = V(1,:)./v;
voy = V(2,:)./v;
voz = V(3,:)./v;

% First we find the Z position where V hits the wall at Rc
theta_f = zeros(1,np);
% discarded = false(1,np);
X_f = zeros(1,np);
Y_f = zeros(1,np);
Z_f = zeros(1,np);

a = vox.^2 + voy.^2;
b = 2*(xo.*vox + yo.*voy);
c = xo.^2 + yo.^2 - Rc^2;
ciw = xo.^2 + yo.^2 - camera_params.Riw^2;

discriminant = b.^2 - 4*a.*c;
discriminantiw = b.^2 - 4*a.*ciw;

INDEX = (discriminant < 0) | (discriminantiw >= 0); % this particles are discarded
II = find(INDEX == false);

xp = 0.5*(-b(II) + sqrt(discriminant(II)))./a(II);
xn = 0.5*(-b(II) - sqrt(discriminant(II)))./a(II);

t = max([xp;xn],[],1);

X_f(II) = xo(II) + vox(II).*t;
Y_f(II) = yo(II) + voy(II).*t;
Z_f(II) = zo(II) + voz(II).*t;

theta_f(II) = atan2(Y_f(II),X_f(II));
theta_f(theta_f < 0) = theta_f(theta_f < 0) + 2*pi;

% Then, we calculate the angle between V and the position of the camera
xc = Rc*cos(theta_f(II));
yc = Rc*sin(theta_f(II));

ax = xc - xo(II);
ay = yc - yo(II);
az = Zc - zo(II);

a = sqrt(ax.^2 + ay.^2 + az.^2);

ax = ax./a;
ay = ay./a;
az = az./a;

psi = acos(ax.*vox(II) + ay.*voy(II) + az.*voz(II))';
visible = psi <= angle(II); % only discarded == false

ip = false(1,np);
ip(II(visible)) = true;

% ip(II) = true;

% Optional figures
if ( option )
    XP = X(:,ip);
    XC = [X_f(ip);Y_f(ip);Z_f(ip)];
    theta = theta_f(ip);
    clockwise_rotation = @(t,x) [cos(t),sin(t);-sin(t),cos(t)]*x;
    for jj=1:size(XP,2)
        XP(1:2,jj) = clockwise_rotation(theta(jj),XP(1:2,jj));
        XC(1:2,jj) = clockwise_rotation(theta(jj),XC(1:2,jj));
    end
    
    Ro = sqrt(xo(ip).^2 + yo(ip).^2);
    figure
    subplot(2,2,1)
    plot(Ro,zo(ip),'g.','MarkerSize',5)
    xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    subplot(2,2,2)
    histogram(zo(ip))
    xlabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$f(Z)$','Interpreter','latex','FontSize',16)
    subplot(2,2,3)
    plot3([XP(1,:);XC(1,:)],[XP(2,:);XC(2,:)],[XP(3,:);XC(3,:)])
    xlabel('$X$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Y$ (m)','Interpreter','latex','FontSize',16)
    axis equal; box on
    subplot(2,2,4)
    plot(theta_f(ip),'b','MarkerSize',5)
    xlabel('Particle number','Interpreter','latex','FontSize',16)
    ylabel('$\theta$ (rad)','Interpreter','latex','FontSize',16)
end
end

function pixel_grid = setupCameraPixelGrid(camera_params,option)
% Here we set-up the pixel grid of the camera in a coordinate system such
% that the x-axis corresponds to the horizontal dimension, and the y-axis
% corresponds to the vertical dimension, with respect to the floor level.
disp(['Setting up camera pixel grid...'])

xmin = -camera_params.size(1)/2;
xmax = camera_params.size(1)/2;
DX = camera_params.size(1)/camera_params.NX;
pixels_x = (xmin+0.5*DX):DX:(xmax-0.5*DX);

ymin = camera_params.position(2) - camera_params.size(2)/2;
ymax = camera_params.position(2) + camera_params.size(2)/2;
DY = camera_params.size(2)/camera_params.NY;
pixels_y = (ymin+0.5*DY):DY:(ymax-0.5*DY);

pixel_grid = struct;
[pixel_grid.X, pixel_grid.Y] = meshgrid(pixels_x,pixels_y);

pixel_grid.xmin = xmin;
pixel_grid.xmax = xmax;

pixel_grid.ymin = ymin;
pixel_grid.ymax = ymax;

pixel_grid.xnodes = pixels_x;
pixel_grid.ynodes = pixels_y;

pixel_grid.xedges = linspace(xmin,xmax,camera_params.NX+1);
pixel_grid.yedges = linspace(ymin,ymax,camera_params.NY+1);

if option
    [X,Y] = ...
        meshgrid(pixel_grid.xedges,pixel_grid.yedges);
    
    h = figure;
    subplot(1,2,1)
    plot(pixel_grid.X,pixel_grid.Y,'rs','MarkerSize',2)
    hold on
    plot(X,Y,'k',X',Y','k')
    hold off
    axis equal
    title('Camera pixel array','Interpreter','latex','FontSize',12)
    xlabel('$X$ (m)','Interpreter','latex','FontSize',14)
    ylabel('$Y$ (m)','Interpreter','latex','FontSize',14)
    
    
    % xy-plane with geometry of camera
    Rmax = ceil(camera_params.position(1));
    R = 0.5:0.5:Rmax;
    t = linspace(0,2*pi,50);
    for ii=1:numel(R)
        x = R(ii)*cos(t);
        y = R(ii)*sin(t);
        
        figure(h)
        subplot(1,2,2)
        hold on
        plot(x,y,'--','Color',[0.4,0.4,0.4])
        hold off
        axis equal
        box on
    end
    
    figure(h)
    subplot(1,2,2)
    hold on
    plot(camera_params.position(1),0,'ko','MarkerSize',5)
    hold off
    
    N = 100;
    ccd_x = linspace(xmin,xmax,N);
    ccd_y = -camera_params.focal_length*ones(1,N);
    
    a = camera_params.incline;
    
    ccd_x_r = ccd_x*cos(a) - ccd_y*sin(a);
    ccd_y_r = ccd_x*sin(a) + ccd_y*cos(a);
    
    ccd_x_r = ccd_x_r + camera_params.position(1);
    
    figure(h)
    subplot(1,2,2)
    hold on
    plot(ccd_x_r,ccd_y_r,'k','LineWidth',2)
    hold off
    
    % Main line of sight
    max_lenght = 2*Rmax;
    
    line_x = zeros(1,N);
    line_y = linspace(-camera_params.focal_length,max_lenght,N);
    
    line_x_r = line_x*cos(a) - line_y*sin(a);
    line_y_r = line_x*sin(a) + line_y*cos(a);
    
    line_x_r = line_x_r + camera_params.position(1);
    
    figure(h)
    subplot(1,2,2)
    hold on
    plot(line_x_r,line_y_r,'r-.','LineWidth',1)
    hold off
    
    % Cone of sight min
    line_x = zeros(1,N);
    line_y = linspace(0,max_lenght+camera_params.focal_length,N);
    
    b = camera_params.horizontal_angle_view;
    
    line_x_r = line_x*cos(b) + line_y*sin(b);
    line_y_r = -line_x*sin(b) + line_y*cos(b);
        
    line_x_r = line_x_r + xmin;
    line_y_r = line_y_r - camera_params.focal_length;
    
    line_x = line_x_r;
    line_y = line_y_r;
    line_x_r = line_x*cos(a) - line_y*sin(a);
    line_y_r = line_x*sin(a) + line_y*cos(a);
    
    line_x_r = line_x_r + camera_params.position(1);
    
    figure(h)
    subplot(1,2,2)
    hold on
    plot(line_x_r,line_y_r,'r','LineWidth',1)
    hold off
    
    % Cone of sight max
    line_x = zeros(1,N);
    line_y = linspace(0,max_lenght+camera_params.focal_length,N);
    
    line_x_r = line_x*cos(b) - line_y*sin(b);
    line_y_r = line_x*sin(b) + line_y*cos(b);
        
    line_x_r = line_x_r + xmax;
    line_y_r = line_y_r - camera_params.focal_length;
    
    line_x = line_x_r;
    line_y = line_y_r;
    line_x_r = line_x*cos(a) - line_y*sin(a);
    line_y_r = line_x*sin(a) + line_y*cos(a);
    
    line_x_r = line_x_r + camera_params.position(1);
    
    figure(h)
    subplot(1,2,2)
    hold on
    plot(line_x_r,line_y_r,'r','LineWidth',1)
    hold off
end


end

function [rotation_angle,ip_in_pixel] = findRotationAngles(X,camera_params,option)
disp('Calculating rotation angles...')
% Here we compute the rotation angles in the toroidal direction based on
% the pixel array of the camera.
np = size(X,2);
Rcam = camera_params.position(1);
Zcam = camera_params.position(2);
incline = camera_params.incline;

xnodes = camera_params.pixel_grid.xnodes;
xedges = camera_params.pixel_grid.xedges;

ynodes = camera_params.pixel_grid.ynodes;
yedges = camera_params.pixel_grid.yedges;

eta_angle = abs(atan2(xnodes,camera_params.focal_length));
beta_angle = zeros(1,camera_params.NX);
for ii=1:camera_params.NX
    if (xedges(ii) < 0)
        beta_angle(ii) = pi/2 - incline - eta_angle(ii);
    else
        beta_angle(ii) = pi/2 - incline + eta_angle(ii);
    end
end

psi_range = atan2(yedges,camera_params.focal_length);

threshold_angle = atan2(camera_params.Riw,-Rcam);
threshold_radius = sqrt( camera_params.Riw^2 + Rcam^2 );

R = sqrt(sum(X(1:2,:).^2,1));
D = sqrt( (X(1,:) - Rcam).^2 + X(2,:).^2 );
psi = -atan2(X(3,:)-Zcam,D);

rotation_angle = cell(2,camera_params.NX,camera_params.NY);
ip_in_pixel = cell(2,camera_params.NX,camera_params.NY);
for ii=1:camera_params.NX
    a = 1 + tan(beta_angle(ii)).^2;
    b = -2*tan(beta_angle(ii)).^2*Rcam;
    c = (Rcam*tan(beta_angle(ii)))^2 - R.^2;
    
    xp = 0.5*(-b + sqrt(b.^2 - 4*a*c))./a;
    xn = 0.5*(-b - sqrt(b.^2 - 4*a*c))./a;

    ip = find(imag(xp) == 0);
    in = find(imag(xn) == 0);
    
    yp = sqrt(R.^2 - xp.^2);
    yn = sqrt(R.^2 - xn.^2);
    
    if ~isempty(ip)
        xtmp = xp(ip) - Rcam;
        ytmp = yp(ip);
        Rtmp = sqrt(xtmp.^2 + ytmp.^2);
        tmp = atan2(ytmp,xtmp);
        
        I = find(tmp > threshold_angle);
        II = find(Rtmp > threshold_radius);
        ip(intersect(I,II)) = [];
    end
    
    if ~isempty(in)
        xtmp = xn(in) - Rcam;
        ytmp = yn(in);
        Rtmp = sqrt(xtmp.^2 + ytmp.^2);
        tmp = atan2(ytmp,xtmp);
        
        I = find(tmp > threshold_angle);
        II = find(Rtmp > threshold_radius);
        in(intersect(I,II)) = [];
    end
    
    if option
        hold on;plot(xp(ip),yp(ip),'k.','MarkerSize',1);hold off
        hold on;plot(xn(in),yn(in),'k.','MarkerSize',1);hold off
%         INDEX = 1;
%         if ~isempty(intersect(ip,INDEX))
%             hold on;plot(xp(INDEX),yp(INDEX),'k.','MarkerSize',8);hold off
%         end
%         if ~isempty(intersect(in,INDEX))
%             hold on;plot(xn(INDEX),yn(INDEX),'k.','MarkerSize',8);hold off
%         end
%         INDEX = 107;
%         if ~isempty(intersect(ip,INDEX))
%             hold on;plot(xp(INDEX),yp(INDEX),'b.','MarkerSize',8);hold off
%         end
%         if ~isempty(intersect(in,INDEX))
%             hold on;plot(xn(INDEX),yn(INDEX),'b.','MarkerSize',8);hold off
%         end
    end
    
    for jj=1:camera_params.NY
        ip_in_pixel_y = find((psi>=psi_range(jj)) & (psi<psi_range(jj+1)));
        
        ip_in_pixel{1,ii,jj} = intersect(ip,ip_in_pixel_y);
        ip_in_pixel{2,ii,jj} = intersect(in,ip_in_pixel_y);
        
        tmp = atan2(yp(ip_in_pixel{1,ii,jj}),xp(ip_in_pixel{1,ii,jj}));
        tmp(tmp<0) = tmp(tmp<0) + 2*pi;
        rotation_angle{1,ii,jj} = tmp;
        
        tmp = atan2(yn(ip_in_pixel{2,ii,jj}),xn(ip_in_pixel{2,ii,jj}));
        tmp(tmp<0) = tmp(tmp<0) + 2*pi;
        rotation_angle{2,ii,jj} = tmp;
    end
    
end

end

function SD = syntheticDiagnosticSynchrotron(ST,compare)
% Synthetic diagnostic of camera for synchrotron radiation.
SD = struct;
disp('Starting synchrotron synthetic diagnostic...')

lch = 1E2;
lchnm = 1E7;
Pch_tot = 1E-9;
Pch_lambda = 1E-5;

% Simulation parameters
num_species = double(ST.params.simulation.num_species);

% Resolution in wavelength
N = 50;
% Wavelength range (450 nm - 950 nm for visible light).
lambda_min = 400E-9;% in meters
lambda_max = 900E-9;% in meters
Dl = (lambda_max - lambda_min)/N;
lambda = zeros(1,N);
for ll=1:N
    lambda(ll) = lambda_min + Dl*(ll-1);
end
% lambda = linspace(lambda_min,lambda_max,N);
lambda = lch*lambda; % in cm

% Camera parameters
camera_params = struct;
camera_params.Riw = 1.0;% inner wall radius in meters
camera_params.NX = 20;
camera_params.NY = 20;
camera_params.size = [0.4,0.4]; % [horizontal size, vertical size] in meters
camera_params.focal_length = 0.5; % In meters
camera_params.position = [2.4,0.0]; % [R,Z] in meters
% The angle defined by the detector plane (pixel array) and the x-axis of a
% coordinate system where phi = 0, the toroidal angle, corresponds to the
% y-axis, and phi = 90 corresponds to the x-axis
camera_params.incline = 55; % in degrees
camera_params.incline = deg2rad(camera_params.incline);
camera_params.horizontal_angle_view = ...
    atan2(0.5*camera_params.size(1),camera_params.focal_length); % in radians
camera_params.vertical_angle_view = ...
    atan2(0.5*camera_params.size(2),camera_params.focal_length); % in radians
camera_params.pixel_grid = setupCameraPixelGrid(camera_params,true);

clockwise_rotation = @(t,x) [cos(t),sin(t);-sin(t),cos(t)]*x;
anticlockwise_rotation = @(t,x) [cos(t),-sin(t);sin(t),cos(t)]*x;

for ss=1:num_species
    q = abs(ST.params.species.q(ss));
    m = ST.params.species.m(ss);
    Ro = ST.params.fields.Ro;
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    bool = pin;% & passing;
    
    X = [];
    V = [];
    gammae = [];
    Prad = [];
    eta = [];
    B = [];
    E = [];
    for ii=1:ST.num_snapshots
        X = [X,ST.data.(['sp' num2str(ss)]).X(:,bool,ii)];
        V = [V,ST.data.(['sp' num2str(ss)]).V(:,bool,ii)];
        gammae = [gammae;ST.data.(['sp' num2str(ss)]).g(bool,ii)];
        eta = [eta;pi*ST.data.(['sp' num2str(ss)]).eta(bool,ii)/180];
        Prad = [Prad;abs(ST.data.(['sp' num2str(ss)]).Prad(bool,ii))];
        B = [B,ST.data.(['sp' num2str(ss)]).B(:,bool,ii)];
        try
            E = [E,ST.data.(['sp' num2str(ss)]).B(:,bool,ii)];
        catch
            E = zeros(size(B));
        end
    end
    numPart = numel(eta);
    
    vec = zeros(numPart,1);
    binormal = zeros(size(B));
    for ii=1:numPart
        VxE = cross(squeeze(V(:,ii)),squeeze(E(:,ii)));
        VxB = cross(squeeze(V(:,ii)),squeeze(B(:,ii)));
        VxVxB = cross(squeeze(V(:,ii)),VxB);
        aux = VxE + VxVxB;
        vec(ii) = sqrt( aux'*aux );
        binormal(:,ii) = VxE + VxVxB;
        binormal(:,ii) = binormal(:,ii)/sqrt(binormal(:,ii)'*binormal(:,ii));
    end
    
    v = squeeze( sqrt( sum(V.^2,1) ) )';
    k = q*vec./(m*gammae.*v.^3); % Actual curvature
    k = k/1E2; % In cm^-1
    
    clear VxE VxB VxVxB aux E B v
    
    % This angle is used as the primer criterium for deciding whether the
    % particle is seen by the camera or not.
    threshold_angle = (3*k*lambda(end)/(2*pi)).^(1/3);
    
    [ip,theta_f] = findVisibleParticles(X,V,threshold_angle,camera_params,false);
    
%     clear threshold_angle
    
    % Here we drop the particles that are not seen by the camera
    theta_f(~ip) = [];
    X(:,~ip) = [];
    V(:,~ip) = [];
    binormal(:,~ip) = [];
    gammae(~ip) = [];
    k(~ip) = [];
    eta(~ip) = [];
    Prad(~ip) = [];
    
    numPart = numel(gammae); % Here we re-calculate the number of visible particles
    
    for ii=1:numPart
        X(1:2,ii) = clockwise_rotation(theta_f(ii),X(1:2,ii));
        V(1:2,ii) = clockwise_rotation(theta_f(ii),V(1:2,ii));
        binormal(1:2,ii) = ...
            clockwise_rotation(theta_f(ii),binormal(1:2,ii));
    end
    
%     clear ip theta_f
    
    [rotation_angles,ip_in_pixel] = findRotationAngles(X,camera_params,true);
    
    % % % % Beyond this point all variables are in cgs units % % % %
    c = 1E2*ST.params.scales.v;
    qe = 3E9*q;
    m = 1E3*m;
    
    X = 1E2*X;
    V = 1E2*V;
%     k = k/lch;
    
    xcam = lch*[camera_params.position(1);0;camera_params.position(2)];
    
    % functions for calculating the radiation power density
    K13 = @(x) besselk(1/3,x);
    K23 = @(x) besselk(2/3,x);
    
    zeta = @(g,p,k,l) (2*pi./(3*l*k.*g.^3)).*(1 + (g.*p).^2).^1.5;
    
    Po = @(g,k,l) 8*pi*c*qe^2./(3*k.*(l*g).^4);

    fx = @(g,p,x) g.*x./sqrt( 1 + (g.*p).^2 );
    P0 = @(g,p,k,l) -(4*pi*c*qe^2./(k.*(l*g).^4)).*(1 + (g.*p).^2).^2;
    arg = @(g,p,k,l,x) 1.5*zeta(g,p,k,l).*(fx(g,p,x) + fx(g,p,x).^3/3);
    P1 = @(g,p,k,l,x) (g.*p).^2.*K13(zeta(g,p,k,l)).*cos( arg(g,p,k,l,x) )./(1 + (g.*p).^2);
    P2 = @(g,p,k,l,x) -0.5*K13(zeta(g,p,k,l)).*(1 + fx(g,p,x).^2).*cos( arg(g,p,k,l,x) );
    P3 = @(g,p,k,l,x) fx(g,p,x).*K23(zeta(g,p,k,l)).*sin( arg(g,p,k,l,x) );
    
    Psyn = ...
        @(g,p,k,l,x) P0(g,p,k,l).*( P1(g,p,k,l,x) + P2(g,p,k,l,x) + P3(g,p,k,l,x));
    
    counter_psi_chi = zeros(camera_params.NX,camera_params.NY);
    counter_psi = zeros(camera_params.NX,camera_params.NY);    
    Ptot_psi_chi = zeros(camera_params.NX,camera_params.NY);    
    P_psi_chi = cell(camera_params.NX,camera_params.NY);
    P_psi = cell(camera_params.NX,camera_params.NY);
    Ptot_psi = zeros(camera_params.NX,camera_params.NY);
    for ii=1:camera_params.NX
        for jj=1:camera_params.NY
            disp(['Pixel (' num2str(ii) ',' num2str(jj) ') out of ' num2str(camera_params.NX*camera_params.NY)])
            P_psi_chi{ii,jj} = zeros(1,N);
            P_psi{ii,jj} = zeros(1,N);
            
            I = [ip_in_pixel{1,ii,jj},ip_in_pixel{2,ii,jj}];
            
            if ~isempty(I)
                % Rotation angles
                anticlockwise = ...
                    [rotation_angles{1,ii,jj},rotation_angles{2,ii,jj}];
                
                clockwise = atan2(X(2,I),X(1,I));
                clockwise(clockwise<0) = clockwise(clockwise<0) + 2*pi;
                
                angle = anticlockwise - clockwise;
                % Rotation angles
                X_pix = zeros(3,numel(I));
                v_unit = zeros(3,numel(I));
                binormal_pix = zeros(3,numel(I));
                gamma_pix = zeros(1,numel(I));
                kappa_pix = zeros(1,numel(I));
                for kk=1:numel(I)
                    X_pix(:,kk) = X(:,I(kk));
                    v_unit(:,kk) = normc(V(:,I(kk)));
                    binormal_pix(:,kk) = binormal(:,I(kk));
                    gamma_pix(kk) = gammae(I(kk));
                    kappa_pix(kk) = k(I(kk));
                end
                xcam_tmp = [xcam(1).*cos(angle) + xcam(2).*sin(angle);...
                    -xcam(1).*sin(angle) + xcam(2).*cos(angle);...
                    xcam(3)*ones(size(angle))];
                n = normc(xcam_tmp - X_pix);
                
                ndotBinormal = dot(n,binormal_pix,1);
                aa = acos(ndotBinormal);
                II = aa > pi/2;
                
                psi = zeros(size(aa));
                psi(II) = aa(II) - pi/2;
                psi(~II) = pi/2 - aa(~II);
                psi = psi;
                
                nperp = normc(n - repmat(ndotBinormal,3,1).*binormal_pix);
                chi = abs( acos(dot(nperp,v_unit,1)) );
                chi = chi;
                                
                for ll=1:N
                    xi = 2*pi./(3*lambda(ll)*kappa_pix.*gamma_pix.^3);
                    D = ( 0.5*(sqrt(4 + (pi./xi).^2) - pi./xi) ).^(1/3);
                    chic = (1./D - D)./gamma_pix;
                    psic = (1.5*kappa_pix*lambda(ll)/pi).^(1/3);
                    

                    
                    reject = (chi > chic) | (psi > psic);
                    II = find(reject == false);
                    
%                     if ~isempty(II)
%                         for kk=1:numel(II)
%                             disp(['(X,Y)=' num2str(ii) ',' num2str(jj) ' ll: ' num2str(ll) ' PSI: ' ...
%                                 num2str(psic(II)) ' CHI: ' num2str(chic(II))...
%                                 ' Psyn: ' num2str( Psyn(gamma_pix(II),psi(II),kappa_pix(II),lambda(ll),chi(II)) )])
%                         end
%                         disp('Stop!')
%                     end
                    
                    Psyn_tmp = Psyn(gamma_pix(II),psi(II),kappa_pix(II),lambda(ll),chi(II));
                    P_psi_chi{ii,jj}(ll) = sum(Psyn_tmp(Psyn_tmp > 0));
                    counter_psi_chi(ii,jj) = counter_psi_chi(ii,jj) + numel(find(Psyn_tmp > 0));
                    
                    Psyn_tmp = ...
                        Po(gamma_pix,kappa_pix,lambda(ll)).*(1+(gamma_pix.*psi).^2).^2.*( K23(zeta(gamma_pix,psi,kappa_pix,lambda(ll))).^2 + ...
                        K13(zeta(gamma_pix,psi,kappa_pix,lambda(ll))).^2.*(gamma_pix.*psi).^2./(1+(gamma_pix.*psi).^2) );
                    P_psi{ii,jj}(ll) = sum(Psyn_tmp);
                end
                
                Ptot_psi_chi(ii,jj) = trapz(lambda,P_psi_chi{ii,jj});
                Ptot_psi(ii,jj) = trapz(lambda,P_psi{ii,jj});
                counter_psi(ii,jj) = numel(I);
            end
        end
    end
    
    
    fh = figure;
    subplot(4,4,[1,2,5,6])
    surfc(camera_params.pixel_grid.ynodes,camera_params.pixel_grid.xnodes,Pch_tot*Ptot_psi_chi,...
        'LineStyle','none')
    colormap(jet(1024))
    box on; axis square; view([90,90])
    ylabel('$x$-axis of detector','FontSize',14,'Interpreter','latex')
    xlabel('$y$-axis of detector','FontSize',14,'Interpreter','latex')
    
    subplot(4,4,[3,4,7,8])
    surfc(camera_params.pixel_grid.ynodes,camera_params.pixel_grid.xnodes,Pch_tot*Ptot_psi,...
        'LineStyle','none')
    colormap(jet(1024))
    box on; axis square; view([90,90])
    ylabel('$x$-axis of detector','FontSize',14,'Interpreter','latex')
    xlabel('$y$-axis of detector','FontSize',14,'Interpreter','latex')
    
    subplot(4,4,[9,13])
    surfc(camera_params.pixel_grid.ynodes,camera_params.pixel_grid.xnodes,counter_psi_chi,...
        'LineStyle','none')
    colormap(jet(1024)); colorbar('Location','southoutside')
    box on; axis square; view([90,90])
    ylabel('$x$-axis of detector','FontSize',14,'Interpreter','latex')
    xlabel('$y$-axis of detector','FontSize',14,'Interpreter','latex')
    title('$P_{syn}(\lambda,\psi,\chi)$','FontSize',14,'Interpreter','latex')
    
    subplot(4,4,[10,14])
    surfc(camera_params.pixel_grid.ynodes,camera_params.pixel_grid.xnodes,counter_psi,...
        'LineStyle','none')
    colormap(jet(1024)); colorbar('Location','southoutside')    
    box on; axis square; view([90,90])
    ylabel('$x$-axis of detector','FontSize',14,'Interpreter','latex')
    xlabel('$y$-axis of detector','FontSize',14,'Interpreter','latex')
    title('$P_{syn}(\lambda,\psi)$','FontSize',14,'Interpreter','latex')
    
    Psyn_mean = zeros(1,N);
    for ii=1:camera_params.NX
        for jj=1:camera_params.NY
            Psyn_mean = Psyn_mean + Pch_lambda*P_psi_chi{ii,jj};
            
            figure(fh)
            subplot(4,4,[11,12])
            hold on
            plot(lchnm*lambda,Pch_lambda*P_psi_chi{ii,jj})
            hold off
        end
    end
%     hold on
%     plot(lchnm*lambda,Psyn_mean/N,'k','LineWidth',3)
%     hold off
    box on
    xlabel('$\lambda$ (nm)','FontSize',14,'Interpreter','latex')
    ylabel('$P_{syn}(\lambda,\psi,\chi)$','FontSize',14,'Interpreter','latex')
    
    Psyn_mean = zeros(1,N);
    for ii=1:camera_params.NX
        for jj=1:camera_params.NY
            Psyn_mean = Psyn_mean + Pch_lambda*P_psi{ii,jj};
            
            figure(fh)
            subplot(4,4,[15,16])
            hold on
            plot(lchnm*lambda,Pch_lambda*P_psi{ii,jj})
            hold off
        end
    end
    hold on
    plot(lchnm*lambda,Psyn_mean/N,'k','LineWidth',3)
    hold off
    box on
    xlabel('$\lambda$ (nm)','FontSize',14,'Interpreter','latex')
    ylabel('$P_{syn}(\lambda,\psi)$','FontSize',14,'Interpreter','latex')
    
   saveas(fh,[ST.path 'SyntheticDiagnostic_ss_' num2str(ss)],'fig')
end
end

function radiationPlane(ST)

num_species = double(ST.params.simulation.num_species);

% clockwise_rotation = @(t,x) [cos(t),sin(t);-sin(t),cos(t)]*x;
clockwise_rotation = ...
    @(t,x) [cos(t).*x(1,:) + sin(t).*x(2,:);-sin(t).*x(1,:) + cos(t).*x(2,:)];
% anticlockwise_rotation = @(t,x) [cos(t),-sin(t);sin(t),cos(t)]*x;

h = figure;
for ss=1:num_species
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    bool = pin;% & passing;
    
    X = [];
    V = [];
    for ii=1:ST.num_snapshots
        X = [X,ST.data.(['sp' num2str(ss)]).X(:,bool,ii)];
        V = [V,ST.data.(['sp' num2str(ss)]).V(:,bool,ii)];
    end
    
    phi = atan2(X(2,:),X(1,:));
    phi(phi<0) = 2*pi + phi(phi<0);
    
    X(1:2,:) = clockwise_rotation(phi,X(1:2,:));
    V(1:2,:) = clockwise_rotation(phi,V(1:2,:));
    V = normc(V);
    
    n = [0,-1,0];
    
    eta = zeros(1,numel(phi));
    for ii=1:numel(phi)
        eta(ii) = acosd(dot(V(:,ii),n));
    end
    
    figure(h)
    subplot(num_species,1,ss)
    histogram(eta,'Normalization','probability')
    title(['$\theta_0=$' num2str(ST.params.species.etao(ss)) ...
        '$^\circ$, $\bar{\eta}=$' num2str(mean(eta))...
        '$^\circ$, $\sigma_\eta=$' num2str(std(eta)) '$^\circ$'],...
        'Interpreter','latex')
end

end

function calculateTemperatureComponents(ST)
num_species = double(ST.params.simulation.num_species);
num_snapshots = double(ST.num_snapshots)
c = ST.params.scales.v;

h = figure;
for ss=1:num_species
    m = ST.params.species.m(ss);
    q = abs(ST.params.species.q(ss));    
    
    Tpar = zeros(1,num_snapshots);
    Tperp = zeros(1,num_snapshots);
    for ii=1:num_snapshots
        Vpar = ST.data.(['sp' num2str(ss)]).V(1,:,ii);
        Upar = mean(Vpar);
        
        Vperp = sqrt(sum(ST.data.(['sp' num2str(ss)]).V(2:3,:,ii).^2,1));
        Uperp = mean(Vperp);
        
        Tpar(ii) = m*mean((Vpar - Upar).^2)/q;
        Tperp(ii) = m*mean((Vperp - Uperp).^2)/q;
    end
    
    subplot(num_species,1,ss)
    plot(ST.time,Tpar,'r',ST.time,Tperp,'b')
    xlim([min(ST.time) max(ST.time)])
    xlabel('Time (s)','Interpreter','latex','FontSize',14)
    ylabel('$T_{RE}$ (eV)','Interpreter','latex','FontSize',14)
end

end

function Psyn = singleParticleSpectrum(ST,lambda,g,eta)
Psyn = zeros(size(lambda));

q = abs(ST.params.species.q(1));
m = ST.params.species.m(1);
c = ST.params.scales.v;
v = c*sqrt(1-1/g^2);
ep = 8.854E-12;% Electric permitivity

Bo = abs(ST.params.fields.Bo);

k = q*Bo*sin(eta)/(g*m*v);

% k = 2*q*sin(eta)/(g*m*v);
l = lambda; 
lc = 4*pi/(3*k*g^3);

z = lc./l;

BK53 = @(x) besselk(5/3,x);
IntBKv = @(nu,x) (pi/sqrt(2))*(1 - 0.25*(4*nu^2 -1))*erfc(sqrt(x)) + ...
    0.25*(4*nu^2 - 1)*sqrt(0.5*pi./x).*exp(-x);

for ii=1:numel(z)
    if (z(ii) < 0.5)
        a = (2.16/2^(2/3))*z(ii)^(1/3);
        Psyn(ii) = integral(BK53,z(ii),a) + IntBKv(5/3,a);
    elseif (z(ii) >= 0.5) && (z(ii) < 2.5)
        a = 0.72*(z(ii) + 1);
        Psyn(ii) = integral(BK53,z(ii),a) + IntBKv(5/3,a);
    else
        Psyn(ii) = IntBKv(5/3,z(ii));
    end
end

Psyn = c*q^2*Psyn./(sqrt(3)*ep*g^2*l.^3);
end

function f = efRE(ST,p,eta)
q = abs(ST.params.species.q(1));
m = ST.params.species.m(1);
c = ST.params.scales.v;

E = ST.params.pdf_params.E;
Z = ST.params.pdf_params.Zeff;
k = ST.params.pdf_params.k;
t = ST.params.pdf_params.t;
xo = (m*c^2/q)/1E6;

fG = @(x,k,t) x.^(k-1).*exp(-x/t)/(gamma(k)*t.^k);

Emin = sqrt(ST.params.pdf_params.min_p^2 + 1);
Emax = sqrt(ST.params.pdf_params.max_p^2 + 1);
x = linspace(Emin,Emax,1E3);
fo = trapz(x,fG(x,k,t/xo));

fE = @(x,k,t,xo) fG(x,k,t/xo)/fo;

A = @(E,Z,p) 2*E*p.^2./((Z+1)*sqrt(p.^2 + 1));
ft = @(E,Z,p,t) 0.5*A(E,Z,p).*exp(A(E,Z,p).*cos(t))./sinh(A(E,Z,p));

Eo = sqrt(p.^2 + 1);

f = fE(Eo,k,t,xo).*ft(E,Z,p,eta);
end

function [I,STDI] = MonteCarloIntegral(ST,nsamples,p_lims,eta_lims)
pmin = min(p_lims);
pmax = max(p_lims);
Emax = sqrt(pmax^2 + 1);
Emin = sqrt(pmin^2 + 1);
etamin = min(eta_lims);
etamax = max(eta_lims);

A = (Emax-Emin)*(etamax-etamin);

P = rand(2,nsamples);
P(1,:) = pmin + (pmax - pmin)*P(1,:); % momentum
P(2,:) = etamin + (etamax - etamin)*P(2,:); % pitch angle
P(3,:) = sqrt(P(1,:).^2 + 1); % gamma factor

l = ST.params.pdf_params.lambda;

f = zeros(1,nsamples);

for ii=1:nsamples
    f(ii) = ...
        efRE(ST,P(1,ii),P(2,ii))*sin(P(2,ii))*singleParticleSpectrum(ST,l,P(3,ii),P(2,ii));
end

f1 = mean(f);
f2 = sum(f.^2)/nsamples;

I = A*f1;
STDI = A*sqrt((f2 - f1^2)/nsamples);
end

function subPhaseSpaceAnalysis(ST)
Np = 150;
Neta = 200;
minLevel = -5;
numLevels = 10;

g = @(p) sqrt(p.^2 + 1);

q = abs(ST.params.species.q(1));
m = ST.params.species.m(1);
c = ST.params.scales.v;
l = ST.params.pdf_params.lambda;

try
    pmin = ST.params.pdf_params.min_p;
    pmax = ST.params.pdf_params.max_p;
catch
    pmin = sqrt(ST.params.pdf_params.min_g^2 - 1);
    pmax = sqrt(ST.params.pdf_params.max_g^2 - 1);
end
Emax = sqrt(pmax^2 + 1);
Emin = sqrt(pmin^2 + 1);
etamin = 0;
etamax = pi;

pAxis = linspace(pmin,pmax,Np);
etaAxis = linspace(etamin,etamax,Neta);

E = (g(pAxis)*m*c^2/q)/1E6; % MeV

fRE = zeros(Np,Neta);
Psyn_sp = zeros(Np,Neta);
ddPsyndpchi = zeros(Np,Neta);
ddPsyndpchi = zeros(Np,Neta);

if (strcmp(ST.params.species.energy_distribution,'AVALANCHE'))
    Ebar = ST.params.pdf_params.Epar/ST.params.pdf_params.Ec;
    Zeff = ST.params.pdf_params.Zeff;
    Ehat = (Ebar - 1)/(1 + Zeff);
    Cz = sqrt(3*(Zeff + 5)/pi)*ST.params.pdf_params.Clog;
    f = @(p,x) (Ehat/Cz)*p.*exp( -p.*(x/Cz + 0.5*Ehat*(1 - x.^2)./x) )./x;
    
    for pp=1:Np
        for ee=1:Neta
            fRE(pp,ee) = f(pAxis(pp),etaAxis(ee));
        end
    end
elseif (strcmp(ST.params.species.energy_distribution,'EXPERIMENTAL-GAMMA'))
    for pp=1:Np
        for ee=1:Neta
            fRE(pp,ee) = efRE(ST,pAxis(pp),etaAxis(ee));
        end
    end
end

for ll=1:numel(l)
    for pp=1:Np
        for ee=1:Neta
            Psyn_sp(pp,ee) = singleParticleSpectrum(ST,l,g(pAxis(pp)),etaAxis(ee));
            ddPsyndpchi(pp,ee) = fRE(pp,ee)*Psyn_sp(pp,ee);
        end
    end
end

[maxval,ieta] = max(ddPsyndpchi,[],2);
[~,ip] = max(maxval);
xo = etaAxis(ieta(ip));
yo = pAxis(ip);

% % Sub-rectangles analysis
AA = ddPsyndpchi;
AA = log10(AA/max(max(AA)));

th = figure('Visible','off');
[C,~] = contour(etaAxis,pAxis,AA,[minLevel minLevel]);

ID = C(1,:)<0;
C(:,ID) = [];
x1 = xo - min(C(1,:));
x2 = max(C(1,:)) - xo;
y1 = max(C(2,:)) - yo;
y2 = yo - min(C(2,:));
xratio = x1/x2;
yratio = y1/y2;

[PT,STDPT] = MonteCarloIntegral(ST,5E4,[pmin pmax],[etamin etamax]);
scales = [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1];
P = zeros(1,numel(scales));
STDP = zeros(1,numel(scales));
boundaryBox = zeros(2,5,numel(scales));
for ii=1:numel(scales)
    etamin = xo - x1*scales(ii);
    etamax = xo + x2*scales(ii);
    pmin = yo - y2*scales(ii);
%     pmax = yo + y1*scales(ii);
    pmax = yo + y1;
    
    boundaryBox(1,:,ii) = [etamin etamax etamax etamin etamin];
    boundaryBox(2,:,ii) = [pmin pmin pmax pmax pmin];
   
    [P(ii),STDP(ii)] = MonteCarloIntegral(ST,5E4,[pmin pmax],[etamin etamax]);
end

boundaryBox(1,:,:) = rad2deg(boundaryBox(1,:,:));
boundaryBox(2,:,:) = (sqrt(boundaryBox(2,:,:).^2 + 1)*m*c^2/q)/1E6;

xAxis = rad2deg(etaAxis);
yAxis = E;
lAxis = l/1E-9;


% % Figures
h = figure;

A = fRE;
A = log10(A/max(max(A)));
cmax = max(max(A));
cmin = minLevel;
x = linspace(0,4,numLevels);
levels = cmin + (cmax-cmin)*tanh(x);

figure(h)
subplot(1,3,1)
contourf(xAxis,yAxis,A,levels,'LineStyle','none')
xlim([0 90]); colormap(jet(1024)); colorbar; caxis([cmin cmax])
xlabel('$\theta$ ($^\circ$)','Interpreter','latex')
ylabel('$\mathcal{E}$ (MeV)','Interpreter','latex')

AA = ddPsyndpchi;
AA = log10(AA/max(max(AA)));
cmax = max(max(AA));
cmin = minLevel;%min(min(A(isfinite(A))));
x = linspace(0,4,numLevels);
levels = cmin + (cmax-cmin)*tanh(x);

figure(h)
subplot(1,3,2)
contourf(xAxis,yAxis,AA,levels,'LineStyle','none')
xlim([0 90]); colormap(jet(1024)); colorbar; caxis([cmin cmax])
title(['$\lambda=$ ' num2str(lAxis) ' nm'],'Interpreter','latex')
xlabel('$\theta$ ($^\circ$)','Interpreter','latex')
ylabel('$\mathcal{E}$ (MeV)','Interpreter','latex')


etao = rad2deg(etaAxis(ieta(ip)));
Eo = E(ip);

figure(h)
subplot(1,3,2)
hold on;plot(etao,Eo,'ko','MarkerFaceColor',[1,0,1],'MarkerSize',8)
for ii=1:numel(scales)
    for jj=1:2
        figure(h)
        subplot(1,3,jj)
        hold on
        if P(ii)/PT < 0.9
            plot(squeeze(boundaryBox(1,:,ii)),squeeze(boundaryBox(2,:,ii)),'--','LineWidth',2)
        else
            plot(squeeze(boundaryBox(1,:,ii)),squeeze(boundaryBox(2,:,ii)),'k-','LineWidth',2)
        end
        hold off
    end
end

figure(h)
subplot(1,3,3)
plot(100*scales,100*P/PT)
grid minor
xlabel('Scale','Interpreter','latex')
ylabel('$P_R/\mathcal{P}_R$','Interpreter','latex')
end

function regionsOfSE(ST,chiAxis,pAxis,fnum)
Np = numel(pAxis);
Nchi = numel(chiAxis);
minLevel = -5;
numLevels = 10;

q = abs(ST.params.species.q(1));
m = ST.params.species.m(1);
c = ST.params.scales.v;

try
    l = ST.params.synthetic_camera_params.lambda;
catch
    l = linspace(3E-6,8E-6,20);
end

g = @(p) sqrt(p.^2 + 1);
eta = @(x) acos(x);

fRE = zeros(size(fnum));
Psyn_sp = zeros(numel(l),Np,Nchi);
ddPsyndpchi = zeros(numel(l),Np,Nchi);
ddPsyndpchi_theory = zeros(numel(l),Np,Nchi);

if (strcmp(ST.params.species.energy_distribution,'AVALANCHE'))
    Ebar = ST.params.pdf_params.Epar/ST.params.pdf_params.Ec;
    Zeff = ST.params.pdf_params.Zeff;
    Ehat = (Ebar - 1)/(1 + Zeff);
    Cz = sqrt(3*(Zeff + 5)/pi)*ST.params.pdf_params.Clog;
    f = @(p,x) (Ehat/Cz)*p.*exp( -p.*(x/Cz + 0.5*Ehat*(1 - x.^2)./x) )./x;
    
    for pp=1:Np
        for cc=1:Nchi
            fRE(pp,cc) = f(pAxis(pp),chiAxis(cc));
        end
    end
elseif (strcmp(ST.params.species.energy_distribution,'EXPERIMENTAL-GAMMA'))
    for pp=1:Np
        for cc=1:Nchi
            fRE(pp,cc) = efRE(ST,pAxis(pp),eta(chiAxis(cc)));
        end
    end
end

for ll=1:numel(l)
    for pp=1:Np
        for cc=1:Nchi
            Psyn_sp(ll,pp,cc) = singleParticleSpectrum(ST,l(ll),g(pAxis(pp)),eta(chiAxis(cc)));
            ddPsyndpchi(ll,pp,cc) = fnum(pp,cc)*Psyn_sp(ll,pp,cc);
            ddPsyndpchi_theory(ll,pp,cc) = fRE(pp,cc)*Psyn_sp(ll,pp,cc);
        end
    end
end

E = (g(pAxis)*m*c^2/q)/1E6; % MeV
% xAxis = chiAxis;
xAxis = rad2deg(acos(chiAxis));
% yAxis = pAxis;
yAxis = E;
lAxis = l/1E-9;


I = floor(linspace(1,numel(l),15));
ntiles = ceil(sqrt(numel(I)));

h = figure;
hh = figure;
for sp=1:numel(I)
    A = squeeze(ddPsyndpchi_theory(I(sp),:,:));

    A = log10(A/max(max(A)));
%     A(A~=0) = log10(A(A~=0)/max(max(A)));
%     A(~isfinite(A)) = 0;

    cmax = max(max(A));
    cmin = minLevel;%min(min(A(isfinite(A))));
    
    x = linspace(0,4,numLevels);
    levels = cmin + (cmax-cmin)*tanh(x);
    
    figure(h)
    subplot(ntiles,ntiles,sp)    
    contourf(xAxis,yAxis,A,levels,'LineStyle','none')
    colormap(jet(1024)); colorbar; caxis([cmin cmax])
    title(['$\lambda=$ ' num2str(lAxis(I(sp))) ' nm'],'Interpreter','latex')
    xlabel('$\theta$ ($^\circ$)','Interpreter','latex')
    ylabel('$\mathcal{E}$ (MeV)','Interpreter','latex')
    
    
    AA = squeeze(ddPsyndpchi(I(sp),:,:));

    AA = log10(AA/max(max(AA)));
%     AA(AA~=0) = log10(AA(AA~=0)/max(max(AA)));
%     AA(~isfinite(AA)) = 0;

    cmax = max(max(AA));
    cmin = minLevel;%min(min(A(isfinite(A))));
    
    x = linspace(0,4,numLevels);
    levels = cmin + (cmax-cmin)*tanh(x);
    
    figure(hh)
    subplot(ntiles,ntiles,sp)    
    contourf(xAxis,yAxis,AA,levels,'LineStyle','none')
    colormap(jet(1024)); colorbar; caxis([cmin cmax])
    title(['$\lambda=$ ' num2str(lAxis(I(sp))) ' nm'],'Interpreter','latex')
    xlabel('$\theta$ ($^\circ$)','Interpreter','latex')
    ylabel('$\mathcal{E}$ (MeV)','Interpreter','latex')
end

A = fRE;
% A = A/max(max(A));
% A(A~=0) = log10(A(A~=0)/max(max(A)));
A = log10(A/max(max(A)));
% A(~isfinite(A)) = 0;

AA = fnum;
% AA = AA/max(max(AA));
% AA(AA~=0) = log10(AA(AA~=0)/max(max(AA)));
AA = log10(AA/max(max(AA)));
% AA(~isfinite(AA)) = 0;

cmax = max(max(A));
cmin = minLevel;
% cmin = min(min(A(isfinite(A))));


x = linspace(0,4,numLevels);
levels = cmin + (cmax-cmin)*tanh(x);

figure(h)
subplot(ntiles,ntiles,ntiles^2)
contourf(xAxis,yAxis,A,levels,'LineStyle','none')
colormap(jet(1024)); colorbar; caxis([cmin cmax])
xlabel('$\theta$ ($^\circ$)','Interpreter','latex')
ylabel('$\mathcal{E}$ (MeV)','Interpreter','latex')

cmax = max(max(AA));
cmin = minLevel;
% cmin = min(min(A(isfinite(A))));

x = linspace(0,4,numLevels);
levels = cmin + (cmax-cmin)*tanh(x);

figure(hh)
subplot(ntiles,ntiles,ntiles^2)
contourf(xAxis,yAxis,AA,levels,'LineStyle','none')
colormap(jet(1024)); colorbar; caxis([cmin cmax])
xlabel('$\theta$ ($^\circ$)','Interpreter','latex')
ylabel('$\mathcal{E}$ (MeV)','Interpreter','latex')


hhh = figure;
for sp=1:numel(I)
    figure(hhh)
    subplot(ntiles,ntiles,sp)
    A = squeeze(Psyn_sp(I(sp),:,:));
    contourf(xAxis,yAxis,A,17,'LineStyle','none')
    colormap(jet(1024));
    colorbar
    title(['$\lambda=$ ' num2str(lAxis(I(sp))) ' nm'],'Interpreter','latex')
    xlabel('$\theta$ ($^\circ$)','Interpreter','latex')
    ylabel('$\mathcal{E}$ (MeV)','Interpreter','latex')
end

% saveas(h,[ST.path 'ddPdpdchi_theory'],'fig')
% saveas(hh,[ST.path 'ddPdpdchi_simulation'],'fig')
% saveas(hhh,[ST.path 'Synchrotron_radiation'],'fig')
end

function SE_phaseSpaceAnalisys(ST)
nbins_p = 75;
nbins_chi = 100;

for ss=1:ST.params.simulation.num_species
    q = abs(ST.params.species.q(ss));
    m = ST.params.species.m(ss);
    c = ST.params.scales.v;
    
    try
        pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
        passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
        bool = pin;% & passing;
    catch
        bool = true(size(ST.data.(['sp' num2str(ss)]).g(:,1)));
    end
    
    g = [];
    eta = [];
    for ii=1:ST.num_snapshots
        g = [g;ST.data.(['sp' num2str(ss)]).g(bool,ii)];
        eta = [eta;ST.data.(['sp' num2str(ss)]).eta(bool,ii)];
    end
    
    E = g*m*c^2;
    E = E/(q*1E6);
    p = sqrt(g.^2 - 1);
    chi = cos(deg2rad(eta));
    
    if (strcmp(ST.params.species.energy_distribution(ss),'AVALANCHE') || strcmp(ST.params.species.energy_distribution(ss),'EXPERIMENTAL-GAMMA'))
        pmax = max([ST.params.pdf_params.max_p max(p)]);
        pmin = min([ST.params.pdf_params.min_p min(p)]);
    end
    
    if (strcmp(ST.params.species.energy_distribution(ss),'HOLLMANN'))
        pmin = sqrt(ST.params.pdf_params.min_g^2 - 1);
        pmax = sqrt(ST.params.pdf_params.max_g^2 - 1);
    end

    
    
    Emin = sqrt(pmin.^2 + 1)*m*c^2/(q*1E6);
    Emax = sqrt(pmax.^2 + 1)*m*c^2/(q*1E6);
    Dp = (pmax-pmin)/nbins_p;
    pAxis = pmin + (0:1:nbins_p-1)*Dp + 0.5*Dp;
    EAxis = sqrt(pAxis.^2 + 1)*m*c^2/(q*1E6);
    
    if (strcmp(ST.params.species.energy_distribution(ss),'AVALANCHE') || strcmp(ST.params.species.energy_distribution(ss),'EXPERIMENTAL-GAMMA')...
            || strcmp(ST.params.species.energy_distribution(ss),'HOLLMANN'))
        pitchmax = max([ST.params.pdf_params.max_pitch_angle max(eta)]);
        pitchmin = min([ST.params.pdf_params.min_pitch_angle min(eta)]);
    end
    Dpitch = (pitchmax - pitchmin)/nbins_chi;
    pitchAxis = pitchmin + (0:1:nbins_chi-1)*Dpitch + 0.5*Dpitch;
    
    chiAxis = cos(deg2rad(pitchAxis));
    
    fRE = zeros(nbins_p,nbins_chi);
    ip = floor((p - pmin)/Dp) + 1;
    ip(ip>nbins_p) = nbins_p;
    ichi = floor((eta - pitchmin)/Dpitch) + 1;
    ichi(ichi>nbins_chi) = nbins_chi;
    
    for pp=1:numel(p)
        fRE(ip(pp),ichi(pp)) = fRE(ip(pp),ichi(pp)) + 1;
    end
    
    S = trapz(fliplr(chiAxis),trapz(pAxis,fRE));
    fRE = fRE/S;
    
    xAxis = pitchAxis;
    yAxis = EAxis;
    %     xAxis = cos(deg2rad(pitchAxis));
    %     yAxis = pAxis;
    
    regionsOfSE(ST,chiAxis,pAxis,fRE);
    
    figure
    subplot(2,1,1)
    histogram2(eta,E,[nbins_chi,nbins_p],'FaceColor','flat','Normalization','probability','LineStyle','none');
    view([0 90]);axis([0 max(eta) 0 max(E)]);axis square;
    cmp = colormap(jet(1024));
    xlabel('$\eta$ ($^\circ$)','FontSize',14,'Interpreter','latex')
    ylabel('$\mathcal{E}$ (MeV)','FontSize',14,'Interpreter','latex')
    
    pperp = sin(deg2rad(eta)).*p;
    pparallel = cos(deg2rad(eta)).*p;
    
    subplot(2,1,2)
    histogram2(pparallel,pperp,[nbins_chi,nbins_p],'FaceColor','flat','Normalization','probability','LineStyle','none');
    view([0 90]);axis([0 max(pparallel) 0 max(pperp)]);axis equal;
    cmp = colormap(jet(1024));
    xlabel('$p_\parallel$ ($mc$)','FontSize',14,'Interpreter','latex')
    ylabel('$p_\perp$ ($mc$)','FontSize',14,'Interpreter','latex')
end

end

function plotEnergyPitchanglePDF(ST)
N = 100; % Energy
M = 100; % Pitch-angle

me = ST.params.scales.m;
qe = ST.params.scales.q;
c = ST.params.scales.v;

for ss=1:ST.params.simulation.num_species 
    try
        pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    catch
        pin = true(1,size(ST.data.(['sp' num2str(ss)]).g,1));
    end
    
    Eo = ST.data.(['sp' num2str(ss)]).g(pin,1)*me*c^2/(qe*1E6);
%     Eo = (ST.data.(['sp' num2str(ss)]).g(pin,1)-1)*me*c^2/(qe*1E6);
    pao = ST.data.(['sp' num2str(ss)]).eta(pin,1);
    
    Ef = ST.data.(['sp' num2str(ss)]).g(pin,end)*me*c^2/(qe*1E6);
%     Ef = (ST.data.(['sp' num2str(ss)]).g(pin,end)-1)*me*c^2/(qe*1E6);
    paf = ST.data.(['sp' num2str(ss)]).eta(pin,end);
    
%     E = (ST.data.(['sp' num2str(ss)]).g(pin,:)-1)'*me*c^2/(qe*1E6);
%     pa = ST.data.(['sp' num2str(ss)]).eta(pin,:)';
%     
%     E_edges = linspace(min([min(Eo(pin)) min(Ef)]),...
%         max([max(Eo(pin)) max(Ef)]),N);
%     pa_edges = linspace(min([min(pao(pin)) min(paf)]),...
%         max([max(pao(pin)) max(paf)]),M);
    
%     E_edges = linspace(0,50,N);
%     pa_edges = linspace(0,25,M);
    
    h1=figure('Visible',ST.visible);
    subplot(3,1,1)
    ho = histogram2(pao,Eo,'Normalization','pdf','DisplayStyle','tile',...
    'LineStyle','none');
    colormap(jet(1024)); caxis([min(min(ho.Values)) max(max(ho.Values))])
    colorbar
    box on; grid minor
    xlabel('$\eta$ ($^\circ$)','FontSize',14,'Interpreter','latex')
    ylabel('$\mathcal{E}$ (MeV)','FontSize',14,'Interpreter','latex')
    
    subplot(3,1,2)
    hf = histogram2(paf,Ef,'Normalization','pdf','DisplayStyle','tile',...
    'LineStyle','none');
    colormap(jet(1024)); caxis([min(min(ho.Values)) max(max(ho.Values))])
    colorbar
    box on; grid minor
    xlabel('$\eta$ ($^\circ$)','FontSize',14,'Interpreter','latex')
    ylabel('$\mathcal{E}$ (MeV)','FontSize',14,'Interpreter','latex')
    
    ho.NumBins = hf.NumBins;
    ho.XBinEdges = hf.XBinEdges;
    ho.YBinEdges = hf.YBinEdges;
    
    xAxis = 0.5*(hf.XBinEdges(2:end) + hf.XBinEdges(1:end-1));
    yAxis = 0.5*(hf.YBinEdges(2:end) + hf.YBinEdges(1:end-1));
    
    A = ho.Values - hf.Values;
    
    subplot(3,1,3,'ClippingStyle','rectangle')
    surf(xAxis,yAxis,A','LineStyle','none')
    colormap(jet(1024)); colorbar;view([0 90]);box on;grid minor
    xlabel('$\eta$ ($^\circ$)','FontSize',14,'Interpreter','latex')
    ylabel('$\mathcal{E}$ (MeV)','FontSize',14,'Interpreter','latex')
    
    subplot(3,1,1)
    axis([hf.XBinEdges(1) hf.XBinEdges(end) hf.YBinEdges(1) hf.YBinEdges(end)])
    subplot(3,1,2)
    axis([hf.XBinEdges(1) hf.XBinEdges(end) hf.YBinEdges(1) hf.YBinEdges(end)])
    subplot(3,1,3)
    axis([hf.XBinEdges(1) hf.XBinEdges(end) hf.YBinEdges(1) hf.YBinEdges(end)])
    
%     g = figure;
%     subplot(2,1,1)
%     hold on
%     %     plot(pa,E,'b')
%     plot(pao,Eo,'k.',paf,Ef,'r.','MarkerSize',8)
%     hold off
%     box on
%     grid on
%     xlabel('$\eta$ ($^\circ$)','FontSize',14,'Interpreter','latex')
%     ylabel('$\mathcal{E}$ (MeV)','FontSize',14,'Interpreter','latex')
%     % legend({['$t=$' num2str(ST.time(1)) ' s'],...
%     %     ['$t=$' num2str(ST.time(end)) ' s']},'Interpreter','latex')
%     
%     
%     u = paf-pao;
%     v = Ef-Eo;
%     
%     subplot(2,1,2)
%     quiver(pao,Eo,u,v,'MaxHeadSize',0.02,'AutoScale','off')
%     box on
%     grid on
%     xlabel('$\eta$ ($^\circ$)','FontSize',14,'Interpreter','latex')
%     ylabel('$\mathcal{E}$ (MeV)','FontSize',14,'Interpreter','latex')
%     saveas(g,[ST.path 'EnergyvsPitch' num2str(ss)],'fig')
end
end

function figuresAPS2017(ST)

if isfield(ST.params.fields,'a')
    t = linspace(0,2*pi,200);
    Rs = ST.params.fields.Ro + ST.params.fields.a*cos(t);
    Zs = ST.params.fields.a*sin(t);
end

legends = cell(1,ST.params.simulation.num_species);
pitch_angle = zeros(1,ST.params.simulation.num_species);
confined = zeros(3,ST.params.simulation.num_species);
confinedParticles = zeros(numel(ST.time),ST.params.simulation.num_species);
for ss=1:ST.params.simulation.num_species
    numConfPart = sum(ST.data.(['sp' num2str(ss)]).flag(:,1),1);
    confinedParticles(:,ss) = 100*sum(ST.data.(['sp' num2str(ss)]).flag,1)/numConfPart;
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    trapped = logical( any(ST.data.(['sp' num2str(ss)]).eta > 90,2) );
    pp = pin & passing;
    pt = pin & trapped;
    
    legends{ss} = ['$\eta_0 =$' num2str(ST.params.species.etao(ss)) '$^\circ$'];
    
    confined(1,ss) = 100*sum(pp)/numConfPart;
    confined(2,ss) = 100*sum(pt)/numConfPart;
    confined(3,ss) = confinedParticles(end,ss);
    
    pitch_angle(ss) = ST.params.species.etao(ss);
end


h1 = figure;
h1.Position(3:4) = [1320,830];
C = colormap(h1,jet(1024));
offset = floor(1024/ST.params.simulation.num_species);
colour = C(1:offset:end,:);

figure(h1)
subplot(2,3,1)
plot(ST.time,confinedParticles)
legend(legends,'interpreter','latex','FontSize',12)
box on;axis on;
xlabel('Time $t$ (sec)','Interpreter','latex','FontSize',16)
ylabel('$\%$ of confined RE','Interpreter','latex','FontSize',16)

figure(h1)
subplot(2,3,2)
plot(pitch_angle,confined,'o-')
box on;axis on;grid on;grid minor;
xlim([0 180])
xlabel('$\theta_0$ ($^\circ$)','Interpreter','latex','FontSize',16)
ylabel('Confined RE ($\%$)','Interpreter','latex','FontSize',16)


for ss=1:ST.params.simulation.num_species   
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    bool = pin & passing;
    
    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,bool,1));
    R = sqrt( sum(X(1:2,:).^2,1) );
    Z = X(3,:);

    figure(h1)
    subplot(2,3,4)
    hold on
    plot(R,Z,'.','MarkerSize',8,'MarkerFaceColor',colour(ss,:),'MarkerEdgeColor',colour(ss,:))
    hold off
end

figure(h1)
subplot(2,3,4)
legend(legends,'interpreter','latex','FontSize',12)
if isfield(ST.params.fields,'a')
    hold on
    plot(Rs,Zs,'k','LineWidth',2,'DisplayName','Plasma edge')
    hold off
end
box on;axis on;grid on;grid minor;axis square
xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)

for ss=ST.params.simulation.num_species:-1:1
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    trapped = logical( any(ST.data.(['sp' num2str(ss)]).eta > 90,2) );
    bool = pin & trapped;
    
    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,bool,1));
    R = sqrt( sum(X(1:2,:).^2,1) );
    Z = X(3,:);
    
    figure(h1)
    subplot(2,3,5)
    hold on
    plot(R,Z,'.','MarkerSize',8,'MarkerFaceColor',colour(ss,:),'MarkerEdgeColor',colour(ss,:))
    hold off
end


figure(h1)
subplot(2,3,5)
legend(fliplr(legends),'interpreter','latex','FontSize',12)
if isfield(ST.params.fields,'a')
    hold on
    plot(Rs,Zs,'k','LineWidth',2,'DisplayName','Plasma edge')
    hold off
end
box on;axis on;grid on;grid minor;axis square
xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)


PR = zeros(3,ST.params.simulation.num_species);
for ss=1:ST.params.simulation.num_species     
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    trapped = logical( any(ST.data.(['sp' num2str(ss)]).eta > 90,2) );
    pp = pin & passing;
    pt = pin & trapped;
    try
        PR(1,ss) = sum(abs(ST.data.(['sp' num2str(ss)]).Prad(pp,end)),1);
        PR(2,ss) = sum(abs(ST.data.(['sp' num2str(ss)]).Prad(pt,end)),1);
        PR(3,ss) = PR(1,ss) + PR(2,ss);
    catch
        disp('MESSAGE: No raw radiation data.')
    end
end

figure(h1)
subplot(2,3,3)
plot(pitch_angle,PR,'o-')
box on;axis on;grid on;grid minor;
xlim([0 180])
xlabel('$\theta_0$ ($^\circ$)','Interpreter','latex','FontSize',16)
ylabel('$\mathcal{P}_R$ (Watts)','Interpreter','latex','FontSize',16)

ep0 = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass
Bo = abs(ST.params.fields.Bo);

hh = figure;
hs1 = figure;

PR_ratio = zeros(3,ST.params.simulation.num_species);
etao = zeros(1,ST.params.simulation.num_species);
Po = zeros(1,ST.params.simulation.num_species);
for ss=1:ST.params.simulation.num_species
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    trapped = logical( any(ST.data.(['sp' num2str(ss)]).eta > 90,2) );
    pp = pin & passing;
    pt = pin & trapped;
    
%     B = sqrt(sum(ST.data.(['sp' num2str(ss)]).B(:,pin,end).^2,1))';
    etao(ss) = ST.params.species.etao(ss);
    go = ST.params.species.go(ss);
    vo = c*sqrt(1-1/go^2);
    g = ST.data.(['sp' num2str(ss)]).g(pin,end);
    v = c*sqrt(1-1./g.^2);
    eta_pp = ST.data.(['sp' num2str(ss)]).eta(pp,end);
    eta_pt = ST.data.(['sp' num2str(ss)]).eta(pt,end);
    eta = ST.data.(['sp' num2str(ss)]).eta(pin,end);
    
    ko = qe*Bo*sin(deg2rad(etao(ss)))/(go*vo*me);
    kapp = qe*Bo*sin(deg2rad(eta))./(g.*v*me);
    
    try
        Po(ss) = qe^2*(go*vo)^4*ko^2/(6*pi*ep0*c^3);
        Papp = qe^2*(g.*v).^4.*kapp.^2./(6*pi*ep0*c^3);
        Prad_pp = abs(ST.data.(['sp' num2str(ss)]).Prad(pp,end));
        Prad_pt = abs(ST.data.(['sp' num2str(ss)]).Prad(pt,end));
        Prad = abs(ST.data.(['sp' num2str(ss)]).Prad(pin,end));
        
        
        hs = figure;
        subplot(4,5,[2 3 4 5 7 8 9 10 12 13 14 15])
        plot(eta_pp,Prad_pp,'k.',eta_pt,Prad_pt,'r.',eta,Papp,'g.','MarkerSize',2)
        hold on;plot(etao(ss),Po(ss),'mo','MarkerSize',6,'MarkerFaceColor',[1,0,1],'MarkerEdgeColor',[0,0,0]);hold off
        box on;grid on;grid minor;
        ax = gca;
        ax.XTickLabel = {};
        ax.YTickLabel = {};
        xlim_main = ax.XLim;
        ylim_main = ax.YLim;
        
        subplot(4,5,[1 6 11])
        hold on;
        histogram(Prad,'normalization','count','LineStyle','none')
        histogram(Prad_pp,'normalization','count','LineStyle','none')
        histogram(Prad_pt,'normalization','count','LineStyle','none')
        hold off
        view([-90 90]);box on;grid on;grid minor;
        ax = gca;ax.YTickLabel = {};xlim(ylim_main)
        hold on;plot(Po(ss)*ones(1,10),linspace(0,ax.YLim(2),10),'r','LineWidth',2);hold off
        xlabel('$\mathcal{P}_R$ (Watts)','Interpreter','latex','FontSize',16)
        
        subplot(4,5,[17 18 19 20])
        hold on
        histogram(eta,'normalization','count','LineStyle','none')
        histogram(eta_pp,'normalization','count','LineStyle','none')
        histogram(eta_pt,'normalization','count','LineStyle','none')
        hold off
        view([0 -90]);box on;grid on;grid minor;
        ax = gca;ax.YTickLabel = {};xlim(xlim_main)
        hold on;plot(etao(ss)*ones(1,10),linspace(0,ax.YLim(2),10),'r','LineWidth',2);hold off
        xlabel('$\theta$ ($^\circ$)','Interpreter','latex','FontSize',16)
        saveas(hs,[ST.path 'Scattered_spp_' num2str(ss)],'fig')
        
        figure(hh)
        hold on;histogram(eta_pp,'normalization','pdf','LineStyle','none');hold off
        box on;grid on;grid minor;
        xlabel('$\theta$ ($^\circ$)','Interpreter','latex','FontSize',16)
        ylabel('$f_{RE}(\theta)$','Interpreter','latex','FontSize',16)
        
        figure(hs1)
        subplot(4,5,[2 3 4 5 7 8 9 10 12 13 14 15])
        hold on;plot(eta,Prad,'.',eta,Papp,'k.','MarkerSize',2);hold off
        hold on;plot(etao(ss),Po(ss),'ro','MarkerSize',6,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0]);hold off
        figure(hs1)
        subplot(4,5,[1 6 11])
        hold on;histogram(Prad,'normalization','pdf','LineStyle','none');hold off
        figure(hs1)
        subplot(4,5,[17 18 19 20])
        hold on;histogram(eta,'normalization','pdf','LineStyle','none');hold off
        
        PR_ratio(1,ss) = mean(Prad_pp)/Po(ss);
        if ~isempty(Prad_pt)
            PR_ratio(2,ss) = mean(Prad_pt)/Po(ss);
        else
            PR_ratio(2,ss) = 0;
        end
        PR_ratio(3,ss) = mean(Prad)/Po(ss);
        
        figure(hs1)
        subplot(4,5,[2 3 4 5 7 8 9 10 12 13 14 15])
        ax = gca;
        ax.XTickLabel = {};
        ax.YTickLabel = {};
        xlim_main = ax.XLim;
        ylim_main = ax.YLim;
        box on;grid on;grid minor;
        
        figure(hs1)
        subplot(4,5,[1 6 11])
        view([-90 90]);box on;grid on;grid minor;
        ax = gca;ax.YTickLabel = {};xlim(ylim_main)
        for ss=1:ST.params.simulation.num_species
            hold on;plot(Po(ss)*ones(1,10),linspace(0,ax.YLim(2),10),'r','LineWidth',2);hold off
        end
        xlabel('$\mathcal{P}_R$ (Watts)','Interpreter','latex','FontSize',16)
        
        figure(hs1)
        subplot(4,5,[17 18 19 20])
        view([0 -90]);box on;grid on;grid minor;
        ax = gca;ax.YTickLabel = {};xlim(xlim_main)
        for ss=1:ST.params.simulation.num_species
            hold on;plot(etao(ss)*ones(1,10),linspace(0,ax.YLim(2),10),'r','LineWidth',2);hold off
        end
        xlabel('$\theta$ ($^\circ$)','Interpreter','latex','FontSize',16)
        
        saveas(hs1,[ST.path 'Scattered_all_spp_' num2str(ss)],'fig')
        
        saveas(hs,[ST.path 'Scattered_spp_' num2str(ss)],'fig')
        
        saveas(hh,[ST.path 'histograms'],'fig')
    catch
        disp('MESSAGE: No raw radiation data.')
    end
end

figure(h1)
subplot(2,3,6)
plot(pitch_angle,PR_ratio(3,:),'o-')
box on;axis on;grid on;grid minor;
xlim([0 180])
xlabel('$\theta_0$ ($^\circ$)','Interpreter','latex','FontSize',16)
ylabel('$\langle \mathcal{P}_R \rangle/P_{app}$','Interpreter','latex','FontSize',16)

 saveas(h1,[ST.path 'FigAPS'],'fig')
end

function NIMROD_figure(ST)

hDE = figure;

for ss=1:ST.params.simulation.num_species   
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    trapped = logical( any(ST.data.(['sp' num2str(ss)]).eta > 90,2) );
    
    bool = pin&passing;
    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,bool,1));
    R = sqrt( sum(X(1:2,:).^2,1) );
    Z = X(3,:);
    
    go = ST.data.(['sp' num2str(ss)]).g(bool,1);
    gf = ST.data.(['sp' num2str(ss)]).g(bool,end);
    DE1 = 100*(go-gf)./go;
    
    figure(hDE)
    subplot(3,1,1)
    histogram(DE1,'normalization','count')
    xlabel('$\Delta \mathcal{E}$ ($\%$)','Interpreter','latex','FontSize',16)
    ylabel('$f_{RE}(\Delta \mathcal{E})$','Interpreter','latex','FontSize',16)
    
    t = linspace(0,2*pi,100);
    x = ST.params.species.Ro(ss) + ST.params.species.ro(ss)*cos(t);
    y = ST.params.species.ro(ss)*sin(t);

    I = true(size(DE1));%(DE1>1.5);
    
    h1 = open('poincare_plot_diverted_1100.fig');
    figure(h1)
    hold on;scatter(R(I),Z(I),12,DE1(I),'filled');hold off
    hold on;plot(x,y,'m','LineWidth',2);hold off
    colormap(jet(1024));colorbar;
    if ~isempty(DE1(I))
        caxis([0 max(DE1(I))]);
    end
    box on;grid on;grid minor;
    xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    
    bool = pin&trapped;
    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,bool,1));
    R = sqrt( sum(X(1:2,:).^2,1) );
    Z = X(3,:);
    
    go = ST.data.(['sp' num2str(ss)]).g(bool,1);
    gf = ST.data.(['sp' num2str(ss)]).g(bool,end);
    DE2 = 100*(go-gf)./go;
    
    figure(hDE)
    subplot(3,1,2)
    histogram(DE2,'normalization','count')
    xlabel('$\Delta \mathcal{E}$ ($\%$)','Interpreter','latex','FontSize',16)
    ylabel('$f_{RE}(\Delta \mathcal{E})$','Interpreter','latex','FontSize',16)
    
    t = linspace(0,2*pi,100);
    x = ST.params.species.Ro(ss) + ST.params.species.ro(ss)*cos(t);
    y = ST.params.species.ro(ss)*sin(t);

    h2 = open('poincare_plot_diverted_1100.fig');
    figure(h2)
    hold on;scatter(R,Z,12,DE2,'filled');hold off
    hold on;plot(x,y,'m','LineWidth',2);hold off
    colormap(jet(1024));colorbar;
    if ~isempty(DE2)
        caxis([0 max(DE2)]);
    end
    box on;grid on;grid minor;
    xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    
    bool = pin;
    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,bool,1));
    R = sqrt( sum(X(1:2,:).^2,1) );
    Z = X(3,:);
    
    go = ST.data.(['sp' num2str(ss)]).g(bool,1);
    gf = ST.data.(['sp' num2str(ss)]).g(bool,end);
    DE3 = 100*(go-gf)./go;
    
    figure(hDE)
    subplot(3,1,3)
    histogram(DE3,'normalization','count')
    xlabel('$\Delta \mathcal{E}$ ($\%$)','Interpreter','latex','FontSize',16)
    ylabel('$f_{RE}(\Delta \mathcal{E})$','Interpreter','latex','FontSize',16)
end

end

function M = movieEnergyPitchAngle(ST)

for ss=1:ST.params.simulation.num_species
    m = ST.params.species.m(ss);
    q = abs(ST.params.species.q(ss));
    c = ST.params.scales.v;
    
    fig=figure;
    fig.Position(3) = 1300;
    M(ST.num_snapshots) = struct('cdata',[],'colormap',[]);
    
    go = ST.params.species.go(ss);
    etao = ST.params.species.etao(ss);
    xio = cos(deg2rad(etao));
    Eo = (go-1)*m*c^2/(q*1E6);
    po = sqrt(go^2 - 1);
    pparo = po*cos(deg2rad(etao));
    ppero = po*sin(deg2rad(etao));
    
    g_max = max(max(ST.data.(['sp' num2str(ss)]).g));
    g_min = min(min(ST.data.(['sp' num2str(ss)]).g));
    
    p_max = sqrt(g_max^2 - 1);
    p_min = sqrt(g_min^2 - 1);
    
    E_max = (g_max-1)*m*c^2/(q*1E6);
    E_min = (g_min-1)*m*c^2/(q*1E6);
    
    eta_max = max(max(ST.data.(['sp' num2str(ss)]).eta));
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    trapped = logical( any(ST.data.(['sp' num2str(ss)]).eta > 90,2) );
    
    inAndPassing = pin&passing;
    inAndTrapped = pin&trapped;
    
    for it=1:ST.num_snapshots
        
        C = logical(all(ST.data.(['sp' num2str(ss)]).flag(:,1:it),2));
        D = ~C;
        P = logical( all(ST.data.(['sp' num2str(ss)]).eta(:,1:it) < 90,2) );
        P = P & C;
        T = logical( any(ST.data.(['sp' num2str(ss)]).eta(:,1:it) > 90,2) );
        T = T & C;
              
        p = sqrt(ST.data.(['sp' num2str(ss)]).g(:,it).^2 - 1);
        E = (ST.data.(['sp' num2str(ss)]).g(:,it)-1)*m*c^2/(q*1E6);
        eta = ST.data.(['sp' num2str(ss)]).eta(:,it);
        
        ppar = p.*cos(deg2rad(eta));
        pper = p.*sin(deg2rad(eta));
        
        xi = cos(deg2rad(eta));
        
        NC = 100*sum(C)/ST.params.species.ppp(ss);
        NP = 100*sum(P)/sum(C);
        NT = 100*sum(T)/sum(C);
        
        figure(fig);subplot(1,3,1)
        plot(etao,Eo,'m*',eta(P),E(P),'k.',eta(T),E(T),'r.',eta(D),E(D),'bx')
        title(['$t=$' num2str(ST.time(it)/1E-3) ' ms'],'Interpreter','latex')
        axis([0 eta_max E_min E_max]);grid minor
        xlabel('$\theta$ ($^\circ$)','Interpreter','latex')
        ylabel('$\mathcal{E}$ (MeV)','Interpreter','latex')
        text(0.1*eta_max,E_min + 0.95*(E_max-E_min),...
            ['C: ' num2str(NC) '$\%$ P: ' num2str(NP) '$\%$ T: ' num2str(NT) '$\%$'],...
            'Interpreter','latex','Color','m','FontSize',12)
        
        figure(fig);subplot(1,3,2)
        plot(pparo,ppero,'m*',ppar(P),pper(P),'k.',ppar(T),pper(T),'r.',ppar(D),pper(D),'bx')
        title(['$t=$' num2str(ST.time(it)/1E-3) ' ms'],'Interpreter','latex')
        axis([-p_max p_max 0 p_max]);grid minor
        xlabel('$p_\parallel$ ($m_ec$)','Interpreter','latex')
        ylabel('$p_\perp$ ($m_ec$)','Interpreter','latex')

        figure(fig);subplot(1,3,3)
        plot(po,xio,'m*',p(P),xi(P),'k.',p(T),xi(T),'r.',p(D),xi(D),'bx')
        title(['$t=$' num2str(ST.time(it)/1E-3) ' ms'],'Interpreter','latex')
        axis([p_min p_max -1 1]);grid minor
        xlabel('$p$ ($m_ec$)','Interpreter','latex')
        ylabel('$\xi$','Interpreter','latex')
        
        drawnow
        M(it) = getframe;
    end
end

end

function P = Psp(ST,Bo,lambda,g,eta)
P = zeros(size(eta));

q = abs(ST.params.species.q(1));
m = ST.params.species.m(1);
c = ST.params.scales.v;

v = c*sqrt(1-1/g^2);
ep = 8.854E-12;% Electric permitivity

k = q*Bo*sin(eta)/(g*m*v);

l = lambda; 
lc = 4*pi./(3*k*g^3);

z = lc./l;

BK53 = @(x) besselk(5/3,x);
IntBKv = @(nu,x) (pi/sqrt(2))*(1 - 0.25*(4*nu^2 -1))*erfc(sqrt(x)) + ...
    0.25*(4*nu^2 - 1)*sqrt(0.5*pi./x).*exp(-x);

for ii=1:numel(z)
    if (z(ii) < 0.5)
        a = (2.16/2^(2/3))*z(ii)^(1/3);
        P(ii) = integral(BK53,z(ii),a) + IntBKv(5/3,a);
    elseif (z(ii) >= 0.5) && (z(ii) < 2.5)
        a = 0.72*(z(ii) + 1);
        P(ii) = integral(BK53,z(ii),a) + IntBKv(5/3,a);
    else
        P(ii) = IntBKv(5/3,z(ii));
    end
end

P = c*q^2*P./(sqrt(3)*ep*g^2*l.^3);
end

function pitchAnglePDFSlices(ST)
slices = 7;
N = 100;

hatE = ST.params.pdf_params.E;
Zeff = ST.params.pdf_params.Zeff;

% Here p is normalized by mc and eta is in radians
A = @(E,Z,p) 2*E*p.^2./((Z+1)*sqrt(p.^2 + 1));
ft = @(E,Z,p,t) 0.5*A(E,Z,p).*exp(A(E,Z,p).*cos(t))./sinh(A(E,Z,p));

for ss=1:ST.params.simulation.num_species
    q = abs(ST.params.species.q(ss));
    m = ST.params.species.m(ss);
    c = ST.params.scales.v;
    
    try
        pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
        passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
        bool = pin;% & passing;
    catch
        bool = true(size(ST.data.(['sp' num2str(ss)]).g(:,1)));
    end
    
    g = [];
    eta = [];
    for ii=1:ST.num_snapshots
        g = [g;ST.data.(['sp' num2str(ss)]).g(bool,ii)];
        eta = [eta;ST.data.(['sp' num2str(ss)]).eta(bool,ii)];
    end
    
    E = (g-1)*m*c^2;
    E = E/(q*1E6);
    p = sqrt(g.^2 - 1);
    chi = cos(deg2rad(eta));
    
    if (strcmp(ST.params.species.energy_distribution(ss),'AVALANCHE') || strcmp(ST.params.species.energy_distribution(ss),'EXPERIMENTAL-GAMMA'))
        pmax = max([ST.params.pdf_params.max_p max(p)]);
        pmin = min([ST.params.pdf_params.min_p min(p)]);
    elseif (strcmp(ST.params.species.energy_distribution(ss),'MONOENERGETIC'))
        pmax = max([sqrt(ST.params.species.go(ss)^2 - 1); max(p)]);
        pmin = min([sqrt(ST.params.species.go(ss)^2 - 1); min(p)]);
    end
    
    Emin = 8.0;
    Emax = 15.0;
    %     Emin = sqrt(pmin.^2 + 1)*m*c^2/(q*1E6);
    %     Emax = sqrt(pmax.^2 + 1)*m*c^2/(q*1E6);
    DE = (Emax-Emin)/slices;
    
    EAxis = linspace(Emin,Emax - DE,slices) + 0.5*DE;
    
    if (~strcmp(ST.params.species.energy_distribution(ss),'MONOPITCH'))
        pitchmax = ST.params.pdf_params.max_pitch_angle;
        pitchmin = ST.params.pdf_params.min_pitch_angle;
    end
    
    pitchAxis = linspace(pitchmin,pitchmax,N);
    chiAxis = cos(deg2rad(pitchAxis));
    
    xAxis = pitchAxis;
    yAxis = EAxis;
    
    fig = figure;
    for sl=1:slices
        go = 1E6*q*EAxis(sl)/(m*c^2);
        po = sqrt(go^2 - 1);
        Psyn = Psp(ST,ST.params.pdf_params.Bo,ST.params.pdf_params.lambda,go,deg2rad(pitchAxis));
        f = Psyn.*ft(hatE,Zeff,po,deg2rad(pitchAxis));
        f = f/trapz(pitchAxis,f);
        
        go = 1E6*q*(EAxis(sl)+0.5*DE)/(m*c^2);
        po = sqrt(go^2 - 1);
        Psyn = Psp(ST,ST.params.pdf_params.Bo,ST.params.pdf_params.lambda,go,deg2rad(pitchAxis));
        fu = Psyn.*ft(hatE,Zeff,po,deg2rad(pitchAxis));
        fu = fu/trapz(pitchAxis,fu);
        
        go = 1E6*q*(EAxis(sl)-0.5*DE)/(m*c^2);
        po = sqrt(go^2 - 1);
        Psyn = Psp(ST,ST.params.pdf_params.Bo,ST.params.pdf_params.lambda,go,deg2rad(pitchAxis));
        fl = Psyn.*ft(hatE,Zeff,po,deg2rad(pitchAxis));
        fl = fl/trapz(pitchAxis,fl);
        
        if (~strcmp(ST.params.species.energy_distribution(ss),'MONOENERGETIC'))
            El = EAxis(sl) - 0.5*DE;
            Eu = EAxis(sl) + 0.5*DE;
            I = (E > El) & (E < Eu);
        else
            El = EAxis(sl);
            Eu = EAxis(sl);
            I = true(size(eta));
        end
        
        pitchAxish = linspace(min(pitchAxis),max(pitchAxis),numel(pitchAxis)+1);
        fh = histcounts(eta(I),pitchAxish);
        fh = fh/trapz(pitchAxish(1:N),fh);
        
        pAh = pitchAxish(1:N) + 0.5*mean(diff(pitchAxish));
        
        mean_pitch = trapz(pitchAxis,pitchAxis.*f);
        mean_pitch_h = trapz(pitchAxish(1:N),pAh.*fh);
        
        figure(fig)
        subplot(5,ceil(slices/5),sl)
        plot(pitchAxis,f,'r',pitchAxis,fl,'b-.',pitchAxis,fu,'b--',pAh,fh,'k')
        legend({'$f_o$,','$f_l$','$f_u$','$f_{\mbox{sim}}$'},'Interpreter','latex')
        %         semilogy(pitchAxis,f,'r',pitchAxis(1:N-1),fh,'b')
        title(['$\mathcal{E}\in($' num2str(El) ',' num2str(Eu)...
            ') MeV $\langle \theta \rangle=$' num2str(mean_pitch)...
            '$^\circ$ $\langle \vartheta \rangle=$' num2str(mean_pitch_h) '$^\circ$'],...
            'Interpreter','latex')
        xlabel('$\theta$ ($^\circ$)','Interpreter','latex')
    end
    saveas(fig,[ST.path 'pitch_angle_slices_spp_' num2str(ss) '_range_' num2str(ST.range(2))],'fig')
end

end

function pitchAnglePDFSlicesVsTime(ST)
slices = 7;
N = 100;

Emin = 8.0;
Emax = 15.0;

DE = (Emax-Emin)/slices;
EAxis = linspace(Emin,Emax - DE,slices) + 0.5*DE;

hatE = ST.params.pdf_params.E;
Zeff = ST.params.pdf_params.Zeff;

% Here p is normalized by mc and eta is in radians
A = @(E,Z,p) 2*E*p.^2./((Z+1)*sqrt(p.^2 + 1));
ft = @(E,Z,p,t) 0.5*A(E,Z,p).*exp(A(E,Z,p).*cos(t))./sinh(A(E,Z,p));

for ss=1:ST.params.simulation.num_species
    q = abs(ST.params.species.q(ss));
    m = ST.params.species.m(ss);
    c = ST.params.scales.v;
    
    E_stats = zeros(2,ST.num_snapshots,slices);
    E_error = zeros(ST.num_snapshots,slices);
    
    pa_stats = zeros(2,ST.num_snapshots,slices);
    pa_error = zeros(ST.num_snapshots,slices);
    
    for ii=1:ST.num_snapshots
        try
            pin = logical(ST.data.(['sp' num2str(ss)]).flag(:,ii));
            passing = logical( ST.data.(['sp' num2str(ss)]).eta(:,ii) < 90 );
            trapped = logical( ST.data.(['sp' num2str(ss)]).eta(:,ii) >= 90 );
            bool = pin;% & passing;
        catch
            bool = true(size(ST.data.(['sp' num2str(ss)]).g(:,1)));
        end
        
        g = ST.data.(['sp' num2str(ss)]).g(bool,ii);
        eta = ST.data.(['sp' num2str(ss)]).eta(bool,ii);
        
        E = (g-1)*m*c^2;
        E = E/(q*1E6);
        
        for sl=1:slices
            El = EAxis(sl) - 0.5*DE;
            Eu = EAxis(sl) + 0.5*DE;
            I = (E > El) & (E < Eu);
            
            E_stats(1,ii,sl) = mean(sin(eta(I)).*E(I));
            E_stats(2,ii,sl) = std(sin(eta(I)).*E(I));
            E_error(ii,sl) = E_stats(2,ii,sl)/sqrt(sum(bool));
            
            pa_stats(1,ii,sl) = mean(sin(eta(I)).*eta(I));
            pa_stats(2,ii,sl) = std(sin(eta(I)).*eta(I));
            pa_error(ii,sl) = pa_stats(2,ii,sl)/sqrt(sum(bool));
        end
    end
    
    fig1 = figure;
    fig2 = figure;
    for sl=1:slices
        El = EAxis(sl) - 0.5*DE;
        Eu = EAxis(sl) + 0.5*DE;
        
        xAxis = [ST.time,flip(ST.time)]/1E-3;
        fu = squeeze(E_stats(1,:,sl)) + squeeze(E_error(:,sl))';
        fl = squeeze(E_stats(1,:,sl)) - squeeze(E_error(:,sl))';
        fx = [fu,flip(fl)];
        
        figure(fig1)
        subplot(5,ceil(slices/5),sl)
        fill(xAxis,fx,[0.6,0.6,0.6],'LineStyle','none')
        hold on;
        plot(ST.time/1E-3,squeeze(E_stats(1,:,sl)),'m','LineWidth',2)
        hold off
        title(['$\mathcal{E}\in($' num2str(El) ',' num2str(Eu) ') MeV'],...
            'Interpreter','latex')
        xlabel('$t$ (ms)','Interpreter','latex')
        ylabel('$\langle \mathcal{E} \rangle$ (MeV)','Interpreter','latex')
        
        fu = squeeze(pa_stats(1,:,sl)) + squeeze(pa_error(:,sl))';
        fl = squeeze(pa_stats(1,:,sl)) - squeeze(pa_error(:,sl))';
        fx = [fu,flip(fl)];
        
        figure(fig2)
        subplot(5,ceil(slices/5),sl)
        fill(xAxis,fx,[0.6,0.6,0.6],'LineStyle','none')
        hold on;
        plot(ST.time/1E-3,squeeze(pa_stats(1,:,sl)),'m','LineWidth',2)
        hold off
        title(['$\mathcal{E}\in($' num2str(El) ',' num2str(Eu) ') MeV'],...
            'Interpreter','latex')
        xlabel('$t$ (ms)','Interpreter','latex')
        ylabel('$\langle \theta \rangle$ ($^\circ$)','Interpreter','latex')
    end
    
    saveas(fig1,[ST.path 'energy_slices_vs_time_spp_' num2str(ss) '_range_' num2str(ST.range(2))],'fig')
    saveas(fig2,[ST.path 'pitch_angle_slices_vs_time_spp_' num2str(ss) '_range_' num2str(ST.range(2))],'fig')
end

end

function PAE_PDF(ST)
NE = 100;
NP = 360;

q = abs(ST.params.species.q(1));
m = ST.params.species.m(1);
c = ST.params.scales.v;

Eo = m*c^2;
E = linspace(1.0E6,20.0E6,NE)*q + Eo;
EAxis = (E - Eo)/(q*1E6);
g = E/(m*c^2);
p = sqrt(g.^2 - 1);

pitch = linspace(0,pi,NP);
pitchAxis = pitch*180/pi;

hatE = ST.params.pdf_params.E;
Zeff = ST.params.pdf_params.Zeff;

% Here p is normalized by mc and eta is in radians
A = @(E,Z,p) 2*E*p.^2./((Z+1)*sqrt(p.^2 + 1));
ft = @(E,Z,p,t) 0.5*A(E,Z,p).*exp(A(E,Z,p).*cos(t))./sinh(A(E,Z,p));

mean_pitch = zeros(1,NE);

for ii=1:NE
    Psyn = Psp(ST,ST.params.pdf_params.Bo,ST.params.pdf_params.lambda,g(ii),pitch);
    f = Psyn.*ft(hatE,Zeff,p(ii),pitch);
    PT = trapz(pitch,f.*sin(pitch));
    mean_pitch(ii) = trapz(pitch, pitch.*f.*sin(pitch)/PT);
end

mean_pitch = mean_pitch*180/pi;

figure
plot(EAxis,mean_pitch)


end

function spatialAndVelocityPDF(ST)
load('C-struct.mat')
num_fs = 7;
q = [1 1.5 2 3 4 5 6];

N = 100; % Energy
M = 100; % Pitch-angle

Eslices = 7;
Emin = 8.0;
Emax = 15.0;

DE = (Emax-Emin)/Eslices;
EAxis = linspace(Emin,Emax - DE,Eslices) + 0.5*DE;

% % Parameters

me = ST.params.scales.m;
qe = ST.params.scales.q;
c = ST.params.scales.v;

Ro = C.Ro;
Zo = C.Zo;

for ss=1:ST.params.simulation.num_species 
    try
        pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    catch
        pin = true(1,size(ST.data.(['sp' num2str(ss)]).g,1));
    end
    
    
    X = ST.data.(['sp' num2str(ss)]).X(:,pin,:);
    X = reshape(X,[3 size(X,2)*size(X,3)]);
    
    R = sqrt(sum(X(1:2,:).^2,1));
    Z = X(3,:);
    
    r = sqrt((R-Ro).^2 + (Z-Zo).^2);
    a = atan2(Z-Zo,R-Ro);
    a(a<0) = a(a<0) + 2*pi;
    
%     E = (ST.data.(['sp' num2str(ss)]).g(pin,:)-1)*me*c^2/(qe*1E6);
    E = ST.data.(['sp' num2str(ss)]).g(pin,:)*me*c^2/(qe*1E6);
    E = reshape(E,[size(E,1)*size(E,2) 1]);
    
    pa = ST.data.(['sp' num2str(ss)]).eta(pin,:);
    pa = reshape(pa,[size(pa,1)*size(pa,2) 1]);
       
    
    npart = numel(pa);
    nparto = 0.9626*double(ST.params.species(ss).ppp*ST.params.simulation.nmpi)*size(ST.data.(['sp' num2str(ss)]).X,3);
    nSnapshots = size(ST.data.(['sp' num2str(ss)]).X,3);
    
    
    fig1 = figure;
    fig2 = figure;
    fige = figure;
    for rr=1:num_fs
        E_in = [];
        pa_in = [];
       
        if rr == 1
            rp = interp1(C.(['q_' num2str(rr)]).a,C.(['q_' num2str(rr)]).r,a,'pchip'); 
            for pp=1:npart
                if (r(pp) < rp(pp))
                    E_in = [E_in E(pp)];
                    pa_in = [pa_in pa(pp)];
                end
            end
        else
            rp1 = interp1(C.(['q_' num2str(rr-1)]).a,C.(['q_' num2str(rr-1)]).r,a,'pchip'); 
            rp2 = interp1(C.(['q_' num2str(rr)]).a,C.(['q_' num2str(rr)]).r,a,'pchip'); 
            for pp=1:npart
                if ((r(pp) < rp2(pp)) & (r(pp) > rp1(pp)))
                    E_in = [E_in E(pp)];
                    pa_in = [pa_in pa(pp)];
                end
            end
        end
        
        figure(fig1)
        subplot(num_fs,1,rr)
        ho = histogram2(pa_in,E_in,'FaceColor','flat',...
            'DisplayStyle','tile','ShowEmptyBins','on','LineStyle','none');
        colormap(jet(1024)); caxis([min(min(ho.Values)) max(max(ho.Values))])
        colorbar; box on; grid minor
        xlabel('$\eta$ ($^\circ$)','FontSize',14,'Interpreter','latex')
        ylabel('$\mathcal{E}$ (MeV)','FontSize',14,'Interpreter','latex')
        title(num2str(sum(sum(ho.BinCounts))))
        
        E_in = [];
        pa_in = [];
        
        rp = interp1(C.(['q_' num2str(rr)]).a,C.(['q_' num2str(rr)]).r,a,'pchip');
        for pp=1:npart
            if (r(pp) < rp(pp))
                E_in = [E_in E(pp)];
                pa_in = [pa_in pa(pp)];
            end
        end
        
        
        figure(fig2)
        subplot(num_fs,1,rr)
        ho = histogram2(pa_in,E_in,'FaceColor','flat',...
            'DisplayStyle','tile','ShowEmptyBins','on','LineStyle','none');
        colormap(jet(1024)); caxis([min(min(ho.Values)) max(max(ho.Values))])
        colorbar; box on; grid minor
        xlabel('$\eta$ ($^\circ$)','FontSize',14,'Interpreter','latex')
        ylabel('$\mathcal{E}$ (MeV)','FontSize',14,'Interpreter','latex')
        title(num2str(sum(sum(ho.BinCounts))))
        
        figure(fige)
        hold on
        plot(C.(['q_' num2str(rr)]).fs(1,:),C.(['q_' num2str(rr)]).fs(2,:))
        hold off
        
    end
    
    figure(fige)
    box on;
    grid minor;
    axis equal;
    
    saveas(fig1,[ST.path 'pdfs_rings_spp_' num2str(ss) '_range_' num2str(ST.range(end))],'fig')
    saveas(fig2,[ST.path 'pdfs_disks_spp_' num2str(ss) '_range_' num2str(ST.range(end))],'fig')
   
    
    np = zeros(2,num_fs);
    
    pa = zeros(2,num_fs);
    Upa = zeros(2,num_fs);
    
    np_slices = zeros(Eslices,num_fs);
    pa_slices = zeros(Eslices,num_fs);
    Upa_slices = zeros(Eslices,num_fs);
    
    E = zeros(2,num_fs);
    UE = zeros(2,num_fs);
    
    V = [C.Vp(1), C.Vp(2:end) - C.Vp(1:end-1)];
    
    for ii=1:num_fs
        figure(fig1)
        subplot(num_fs,1,ii)
        ax = gca;
        
        np(1,ii) = sum(sum(ax.Children.BinCounts));
        pa(1,ii) = mean(ax.Children.Data(:,1));
        Upa(1,ii) = std(ax.Children.Data(:,1))/sqrt(np(1,ii));
        E(1,ii) = mean(ax.Children.Data(:,2));
        UE(1,ii) = std(ax.Children.Data(:,2))/sqrt(np(1,ii));
        
        for sl=1:Eslices
            El = EAxis(sl) - 0.5*DE;
            Eu = EAxis(sl) + 0.5*DE;
            I = (ax.Children.Data(:,2) > El) & (ax.Children.Data(:,2) < Eu);
            
            np_slices(sl,ii) = sum(I);
            
            pa_slices(sl,ii) = mean(ax.Children.Data(I,1));
            Upa_slices(sl,ii) = std(ax.Children.Data(I,1))/sqrt(np_slices(sl,ii));
            
            np_slices(sl,ii) = (np_slices(sl,ii)/nSnapshots)./V(ii);
        end
        
        figure(fig2)
        subplot(num_fs,1,ii)
        ax = gca;
        
        np(2,ii) = sum(sum(ax.Children.BinCounts));
        pa(2,ii) = mean(ax.Children.Data(:,1));
        Upa(2,ii) = std(ax.Children.Data(:,1))/sqrt(np(2,ii));
        E(2,ii) = mean(ax.Children.Data(:,2));
        UE(2,ii) = std(ax.Children.Data(:,2))/sqrt(np(2,ii));
    end 
       
    np(1,:) = (np(1,:)/nSnapshots)./V;
    np(2,:) = 100*np(2,:)/nparto;
    
    fig3 = figure;
    subplot(2,3,1)
    plot(q,np(1,:),'o-')
    for sl=1:Eslices
        hold on;
        plot(q,np_slices(sl,:),'o-')
        hold off
    end
    box on;grid minor;
    xlabel('$q$','Interpreter','latex')
    ylabel('$Number of RE$','Interpreter','latex')
    
    subplot(2,3,2)
    errorbar(q,E(1,:),UE(1,:),'o-')
    box on;grid minor;
    xlabel('$q$','Interpreter','latex')
    ylabel('$\left\langle \mathcal{E} \right\rangle$','Interpreter','latex')
    
    subplot(2,3,3)
    errorbar(q,pa(1,:),Upa(1,:),'o-')
    for sl=1:Eslices
        hold on;
        errorbar(q,pa_slices(sl,:),Upa_slices(sl,:),'o-')
        hold off
    end
    box on;grid minor;
    xlabel('$q$','Interpreter','latex')
    ylabel('$\left\langle \theta \right\rangle$','Interpreter','latex')
    
    subplot(2,3,4)
    plot(q,np(2,:),'o-')
    box on;grid minor;
    xlabel('$q$','Interpreter','latex')
    ylabel('$Number of RE$','Interpreter','latex')
    
    subplot(2,3,5)
    errorbar(q,E(2,:),UE(2,:),'o-')
    box on;grid minor;
    xlabel('$q$','Interpreter','latex')
    ylabel('$\left\langle \mathcal{E} \right\rangle$','Interpreter','latex')
    
    subplot(2,3,6)
    errorbar(q,pa(2,:),Upa(2,:),'o-')
    box on;grid minor;
    xlabel('$q$','Interpreter','latex')
    ylabel('$\left\langle \theta \right\rangle$','Interpreter','latex')
    
    saveas(fig3,[ST.path 'pdfs_stats_spp_' num2str(ss) '_range_' num2str(ST.range(end))],'fig')
end

end

function plotEnergyPDF(ST)
N = 100; % Energy bins

me = ST.params.scales.m;
qe = ST.params.scales.q;
c = ST.params.scales.v;

for ss=1:ST.params.simulation.num_species 
    try
        pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    catch
        pin = true(1,size(ST.data.(['sp' num2str(ss)]).g,1));
    end
    
    Eo = ST.data.(['sp' num2str(ss)]).g(pin,1)*me*c^2/(qe*1E6);
%     Eo = (ST.data.(['sp' num2str(ss)]).g(pin,1)-1)*me*c^2/(qe*1E6);
    pao = ST.data.(['sp' num2str(ss)]).eta(pin,1);
    
    Ef = ST.data.(['sp' num2str(ss)]).g(pin,end)*me*c^2/(qe*1E6);
%     Ef = (ST.data.(['sp' num2str(ss)]).g(pin,end)-1)*me*c^2/(qe*1E6);
    paf = ST.data.(['sp' num2str(ss)]).eta(pin,end);
    
    
    fig = figure('Visible','off');
    
    figure(fig)
    subplot(2,1,1)
    histogram(Eo,'Normalization','probability')
    xlabel('Energy (MeV)','Interpreter','latex')
    ylabel('$f_{RE}(\mathcal{E})$','Interpreter','latex')
    
    figure(fig)
    subplot(2,1,2)
    histogram(Ef,'Normalization','probability')
    xlabel('Energy (MeV)','Interpreter','latex')
    ylabel('$f_{RE}(\mathcal{E})$','Interpreter','latex')
end
end
