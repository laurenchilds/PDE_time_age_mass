clear all; close all; clc


%% Data Files

filenames = {...
    'Data/Output_IC2.mat', ... %1
    'Data/Output_IC3.mat', ... %2
    'Data/Output_D1.mat', ... %3 
    'Data/Output_D2.mat', ... %4
    'Data/Output_D3.mat', ... %5
    'Data/Output_D4.mat', ... %6
    'Data/Output_LC1.mat', ... %7
    'placeholder', ... %8
    'placeholder',... %9
    'placeholder', ... %10
    'Data/Output_Asym1.mat',... %11
    'Data/Output_TC.mat',... %12
    'placeholder', ... %13
    'placeholder', ... %14
    'Data/Output_Asym3.mat',... %15
    'Data/Output_Asym2.mat',... %16
    'Data/Output_IC2std-1.mat', ... %17
    'Data/Output_IC2std-2.mat', ... %18
    'Data/Output_D1_c.mat'      ... %19
    'Data/Output_D2_c.mat',     ... %20
    'Data/Output_D3_c.mat',     ... %21
    'Data/Output_D4_c.mat',     ... %22
    'Data/Output_D1_m.mat'      ... %23
    'Data/Output_D2_m.mat',     ... %24
    'Data/Output_D3_m.mat',     ... %25
    'Data/Output_D4_m.mat',         ... %26
    'Data/Output_LC1_v2.mat',... %27
    'Data/Output_LC2_v2.mat',... %28
    'Data/Output_IC3_std_high.mat',... %29
    'Data/Output_IC3_std_low.mat'};     %30
%% Versions
%{
% If not stated, 
%    mu = [0.4, 0.6] 
%    sigma = [0.1 0; 0 0.1] 
%    delta = 0

% Version 1: gTA = 1, gM = m/2, mu = [0.2, 0.6], sigma = [0.1 0; 0 0.1]
% Version 2: gTA = 1, gM = m/2, mu = [0.6, 0.2], sigma = [0.1 0; 0 0.1]
% Version 3: gTA = 1, gM = 1, delta = 0
% Version 4: gTA = 1, gM = 1, delta = a/2
% Version 5: gTA = 1, gM = 1, delta = m/2
% Version 6: gTA = 1, gM = 1, delta = am/4
% Version 7: gTA = 1, gM = m
% Version 8: Does not exist
% Version 9: Does not exist
% Version 10: Does not exit
% Version 11: gTA = 1/10, gM = 1/2(1-m/2)
% Version 12: gTA = 1, gM = m^2/(m^2+4)
% Version 13: Does not exist
% Version 14: Does not exit
% Version 15: gTA = a/10, gM = 1-m^2/(m^2+1)
% Version 16: gTA = 1/10, gM = 1-m^2/(m^2+1)
% Version 17: gTA = 1, gM = m/2, mu = [0.2, 0.6], sigma = [0.1 0; 0 0.5]
% Version 18: gTA = 1, gM = m/2, mu = [0.2, 0.6], sigma = [0.1 0; 0 0.01]
% Version 19: gTA = 1, gM = 1/2, delta = 0
% Version 20: gTA = 1, gM = 1/2, delta = a/2
% Version 21: gTA = 1, gM = 1/2, delta = m/2
% Version 22: gTA = 1, gM = 1/2, delta = am/4
% Version 23: gTA = 1, gM = m/5, delta = 0
% Version 24: gTA = 1, gM = m/5, delta = a/2
% Version 25: gTA = 1, gM = m/5, delta = m/2
% Version 26: gTA = 1, gM = m/5, delta = am/4
% Version 27: gTA = 1/2, gM = m
% Version 28: gTA = a, gM = 1/2
% Version 29: gTA = 1, gM = m/2, mu = [0.6, 0.2], sigma = [0.1 0; 0 0.5]
% Version 30: gTA = 1, gM = m/2, mu = [0.6, 0.2], sigma = [0.1 0; 0 0.01]
%}

%% Versions for Figures

%{
% Figure  1a - Version 3
%         1b - Version 3
%         1c - Version 3
%         1d - Version 3
%         1e - Version 3
% Figure  2a - Version 7
%         2b - Version 12
%         2c - Version 7
%         2d - Version 12        
%         2e - Version 7
%         2f - Version 12 
% Figure  3a - Version 12
%         3b - Version 12
%         3c - Version 12
%         3d - Version 12 
% Figure  4a - Version 27
%         4b - Version 28
%         4c - Version 27
%         4d - Version 28        
%         4e - Version 27
%         4f - Version 28 
% Figure  5a - Version 11
%         5b - Version 16
%         5c - Version 15
%         5d - Version 11        
%         5e - Version 16
%         5f - Version 15
% Figure  6a - Version 3, Version 4, Version 5, Version 6 (top to bottom)
%         6b - Version 19, Version 20, Version 21, Version 22 (top to bottom)
%         6c - Version 23, Version 24, Version 25, Version 26 (top to bottom)
% Figure  7  - Versions 21 and 22
% Figure  8a - Version 18
%         8b - Version 20
%         8c - Version 1
%         8d - Version 2        
%         8e - Version 17
%         8f - Version 29
%}


%% Create Plots

for ver = [ 3 ...               % Figure 1
            7 12 ...            % Figure 2 and Figure 3
            27 28 ...           % Figure 4
            11 16 15 ...        % Figure 5
            3:6 19:26 ...       % Figure 6
            21 22 ...           % Figure 7
            18 30 1 2 17 29]	% Figure 8
    %}
    load(filenames{ver})
    
%%
la = length(avec);
lm = length(mvec);
[A,M] = meshgrid(avec,mvec);

scalef = max(f(:)); % scale to the maximum f (which is density and could be > initial_density)
clear fta ftm ft

fta = zeros(length(tvec),length(avec));
ftm = NaN(length(tvec),length(mvec));
ft = zeros(1,length(tvec));

% Integrate density f(t,a,m) over mass m:
for tcounter = 1:length(tvec)
    for acounter = tcounter:length(avec)
        fta(tcounter,acounter,:) = trapz(mvec,f(tcounter,acounter,:));
    end    
end

% Integrate density f(t,a,m) over age a:
for tcounter = 1:length(tvec)-1
    for mcounter = 1:length(mvec)
        ftm(tcounter,mcounter) = trapz(avec(tcounter:end),f(tcounter,tcounter:end,mcounter));               
    end    
end
ftm(length(tvec),:) = f(length(tvec),end,:)*t_step;

% Integrate density f(t,a,m) over age a and mass m (i.e. integrate fta over a):
for tcounter = 1:length(tvec)
        ft(tcounter) = trapz(avec,fta(tcounter,:));               
end

f1 = f;
for ts = 1:length(tvec)
    for as = 1:length(avec)
        if tvec(ts)>avec(as)
        f1(ts,as,:) = NaN;
        end
    end
end

%% Time Density Plots - Figure 3abcd

if ver == 12

tj = [1 20 40 60];

for j = 1:4
    j1 = tj(j);
figTC = figure(j);


f2 = (squeeze(f1(j1,:,:))');
h = imagesc(avec,mvec,f2);
set(h,'AlphaData',~isnan(f2))
set(gca,'YDir','normal')

xlabel('Age')

if j==1
ylabel('Mass')
end

set(gca,'fontsize',30)
xlim([0 round(a_max)])
ylim([0 round(m_max)])
set(gca,'ytick',[0 round(m_max)],'yticklabel',([0 round(m_max)]))
set(gca,'xtick',[0 round(a_max)],'xticklabel',[0 round(a_max)])
shading interp

saveas(figTC,sprintf('Figures/TimeCourse_%d.fig',j))
saveas(figTC,sprintf('Figures/TimeCourse_%d.png',j))


end
end

%% Total density plot - Figures 1e, 2ef, 4ef

if ver==3 || ver==7 || ver==12 || ver==27 || ver==28
fig1 = figure(100+ver);

% Plot total density (int of f over age and size) as function of time
plot(tvec,movmean(ft,5),'linewidth',2)
ylim([0 1.05])
xlabel('Time, $t$','interpreter','latex')
ylabel('$f(t)$','interpreter','latex')
set(gca,'fontsize',20)
xlim([0 round(t_max)])

saveas(fig1,sprintf('Figures/DensityALL_%d.fig',ver))
saveas(fig1,sprintf('Figures/DensityALL_%d.png',ver))

end

%% Age density histo - Figure 1ac

% Plot total density fta as function of time and age
% Note: each row is a time point, each column is an age group
c = cool(length(tvec));
if ver ==3
fig2 = figure(200+ver);
fta1 = reshape(fta,length(tvec),length(avec));

t1 = find(tvec==1);
t_its = [1 t1:10:length(tvec)];
for j = t_its
h = area(avec,fta1(j,:));
h.FaceColor = c(j,:);
h.FaceAlpha = 0.2;
hold on
tleg{j} = sprintf('t = %1.2f',tvec(j));
end
hold off
xlabel('Age, $a$','interpreter','latex')
ylabel('Density $f(t,a)$','interpreter','latex')
l = legend(tleg{t_its},'fontsize',10);
set(l,'box','off')
set(gca,'fontsize',20)
xlim([0 round(t_max)])

saveas(fig2,sprintf('Figures/HistoAge_%d.fig',ver))
saveas(fig2,sprintf('Figures/HistoAge_%d.png',ver))


end

%% Mass density histo - Figures 1d, 2cd, 4cd, 5def, 7

ftm2 = squeeze(ftm);
if ver == 3 || ver==7 || ver==12 || ver==27 || ver== 28 || ver==16 || ver==11 || ver ==15 || ver==21 || ver==22

% Plot total density ftm as function of time and mass
% Note: each row is a time point, each column is a mass group
fig3 = figure(300+ver);

if ver==21 || ver==22
    t_its = [1 10:10:40];
    fig3 = figure(321);
    hold on
for j = t_its

if ver==21
h = area(mvec,ftm2(j,:));
h.FaceColor = c(1,:);
h.FaceAlpha = 0.1;
else
h = area(mvec,ftm2(j,:));
h.FaceColor = c(end,:);
h.FaceAlpha = 0.1;
end

hold on
tleg{j} = sprintf('t = %1.2f',tvec(j));
end

hold off
xlabel('Mass, $m$','interpreter','latex')
ylabel('Density $f(t,m)$','interpreter','latex')

set(gca,'fontsize',20)
xlim([0 3])
text(0.5,1.37,'t = 0.1')
text(.97,1.17,'t = 1')
text(1.42,.66,'t = 2')
text(2,.25,'t = 3')
text(2.5,.1,'t = 4')

if ver==21
    hleg(1) = h(1);
else
    hleg(2) = h(1);
    l = legend(hleg,{'\delta = m/2','\delta = am/4'});
    set(l,'box','off')
end

else
t1 = find(tvec==1);
t_its = [1 t1:10:length(tvec)];    

for j = t_its
h = area(mvec,ftm2(j,:));
h.FaceColor = c(j,:);
h.FaceAlpha = 0.1;
hold on
tleg{j} = sprintf('t = %1.2f',tvec(j));
end

hold off
xlabel('Mass, $m$','interpreter','latex')
ylabel('Density $f(t,m)$','interpreter','latex')

if ver ~= 3
l = legend(tleg{t_its},'fontsize',10);
set(l,'box','off')

if ver == 11
    set(l,'location','northwest')
    hold on
    a = area([2 2.5],[8 8],'EdgeColor','none');
    a.FaceAlpha = 0.2;
    newcolors = [0.7 0.7 0.7];
    colororder(newcolors)
    a=get(gca,'Children');
    legend(a(2:end),tleg{t_its})
    
end
end
set(gca,'fontsize',20)


xlim([0 round(m_max)])
end
if ver==16 || ver==11 || ver ==15
    xlim([0 2.5])
    set(gca,'xtick',0:2)
end

saveas(fig3,sprintf('Figures/HistoMass_%d.fig',ver))
saveas(fig3,sprintf('Figures/HistoMass_%d.png',ver))


end


%% Heatmap mass - Figures 1b, 2ab, 4ab, 5abc, 6a-c, 8a-f


fig4 = figure(400+ver);

ftm3 = ftm2;

h = imagesc(tvec,mvec,(ftm'));
set(h,'AlphaData',~isnan((ftm')))
set(gca,'fontsize',20)
xlabel('Time, $t$','interpreter','latex')
ylabel('Mass $m$','interpreter','latex')
ylim([0 5])
xlim([0 round(t_max)])
set(gca,'ytick',0:m_max,'yticklabel',(0:m_max))
set(gca,'YDir','normal')
if ver == 3
cbh= colorbar;
cbh.Ticks = cbh.Limits;
cbh.TickLabels = {'Low','High'};
end

if ver==16 || ver==11 || ver ==15
    ylim([0 2.5])
end

saveas(fig4,sprintf('Figures/DensityMass_%d.fig',ver))
saveas(fig4,sprintf('Figures/DensityMass_%d.png',ver))

%% Heatmap Age - Figure 3a

if ver ==3
    
fig3 = figure(300+ver);
fta2 = fta1;

fta3 = fta2;
for ts = 1:length(tvec)
    for as = 1:length(avec)
        if tvec(ts)>avec(as)
        fta3(ts,as) = NaN;
        end
    end
end

h = imagesc(tvec,avec,(fta3'));
set(h,'AlphaData',~isnan((fta3')))
set(gca,'YDir','normal')

set(gca,'fontsize',20)
xlabel('Time, $t$','interpreter','latex')
ylabel('Age $a$','interpreter','latex')
ylim([0 round(t_max)])
xlim([0 round(t_max)])
set(gca,'ytick',0:a_max,'yticklabel',(0:a_max))

end

saveas(fig3,sprintf('Figures/DensityAge_%d.fig',ver))
saveas(fig3,sprintf('Figures/DensityAge_%d.png',ver))

disp(ver)
end



