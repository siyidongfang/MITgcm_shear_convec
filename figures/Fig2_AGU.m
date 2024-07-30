
clear;close all;
addpath ../analysis/colormaps/
fontsize = 19;
load_colors;

addpath ../analysis/
addpath ../analysis/functions/
% expname = 'topo0_H500_s0.0016dz1dx3ln200n-20sm100_kv2e-4';
% expdir = '../exps_hires/';
% expname = 'hires_topo4_s0.0013_dz1dx6n-20';
expname = 'topo4_H500_smo100m_s0.0014_dz1dx3ln200n-20'
expdir = '../backup_exps/exps_topo4/';
% expdir = '/Volumes/MIT/MITgcm_shear_convec/exps_topo4/';
loadexp;
rhoConst = 999.8;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(1)); 
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

%--- snapshots
o = 12*26+3;
tt = squeeze(rdmds([exppath,'/results/THETA'],dumpIters(o)));
tt(tt==0)=NaN;
Hz = sum(delR);
N2const = (1e-3)^2;
tNorth = N2const *(zz+Hz) /9.81/2e-4;
tt_background = ones(Nx,Nr);
for k=1:Nr
    tt_background(:,k) = squeeze(tt_background(:,k))*tNorth(k);
end
tt = tt + tt_background;

rho = rhoConst.*(1-(tt-tRef)*tAlpha);
N2 = NaN*zeros(Nx,Nr);
N2(:,1:Nr-1) = -gravity/rhoConst.*(rho(:,1:end-1)-rho(:,2:end))./(zz(1:end-1)-zz(2:end));
% N2 = N2+N2const;

%--- Load MITgcm simulation
filename = [expdir expname '/RMSE.mat'];
load(filename)
load('fig2/fig2_new.mat')
YLIM = [0 300];

% 
mycolor=cmocean('balance');
mycolor=mycolor(20:end-20,:);
% mycolor=WhiteBlueGreenYellowRed(0);
% mycolor=cmocean('haline');

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1000 300]);

tidx = 120:360;
%%% N2
ax5 = subplot('position',[0.07 0.2 0.85 0.7]);
pcolor(time_tidal(tidx)-time_tidal(tidx(1)),zz-botZ,N2_timeseries(tidx,:)')
hold on;
contour(time_tidal(tidx)-time_tidal(tidx(1)),zz-botZ,N2_timeseries(tidx,:)',[0 0],'Color','c','LineWidth',1);
contour(time_tidal(tidx)-time_tidal(tidx(1)),zz-botZ,uu_timeseries(tidx,:)',[0.15:0.15:0.75],'color',darkgray)
contour(time_tidal(tidx)-time_tidal(tidx(1)),zz-botZ,uu_timeseries(tidx,:)',[0 0],'color',darkgray,'LineWidth',1)
contour(time_tidal(tidx)-time_tidal(tidx(1)),zz-botZ,uu_timeseries(tidx,:)',[-0.75:0.15:-0.15],'--','color',darkgray)
shading interp;
xlabel('Time (tidal cycles)','interpreter','latex');
ylabel('Height above bottom (m)','interpreter','latex');
set(gca,'Fontsize',fontsize);
title('Vertical buoyancy gradient (time series)','Fontsize',fontsize+7,'interpreter','latex','Position',[10,298])
clim(([-1 1]+1)/1e6)
xlim([0 20])
ylim(YLIM)
colormap(mycolor);
h5=colorbar;
set(h5,'Position', [0.94 0.2 0.01 0.700]);
set(get(h5,'Title'),'String',{'$(\mathrm{s}^{-2})\ \ \ \ \  $'},'interpreter','latex','FontSize',fontsize);

print('-dpng','-r300','fig2/AGU1.png');




figure(2)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1000 300]);

ax7 = subplot('position',[0.07 0.2 0.85 0.7]);
pcolor(xx,zz-botZ,N2')
hold on;
contour(xx,zz-botZ,N2',[0 0],'Color','c','LineWidth',1);
hold off;
shading interp;set(gca,'Fontsize',fontsize);
clim(([-1 1]+1)/1e6)
ylabel('Height above bottom (m)','interpreter','latex');
set(gca,'Fontsize',fontsize);
ylim(YLIM)
h7=colorbar(ax7);
set(h7,'Position', [0.94 0.2 0.01 0.700]);
set(get(h7,'Title'),'String',{'$(\mathrm{s}^{-2})\ \ \ \ \  $'},'interpreter','latex','FontSize',fontsize);
xlabel('Across-slope distance (m)','interpreter','latex','FontSize',fontsize+2);
title('Vertical buoyancy gradient (snapshot)','Fontsize',fontsize+7,'interpreter','latex','Position',[0 297])
xlim([-1.5 1.5]*1000)
colormap(mycolor);
print('-dpng','-r300','fig2/AGU2.png');


%%% Save the figure

% print('-dpng','-r300','fig2/fig2.png');