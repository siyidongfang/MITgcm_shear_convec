
clear;close all;
addpath ../analysis/colormaps/
fontsize = 15;
load_colors;

addpath ../analysis/
addpath ../analysis/functions/
expname='topo4_H500Lx3k_s1.70dz1dx3sm100'
expdir = '../exps_kv5e-6/';
topo=4;

% expdir = '../exps_topo4_hires/';
% expdir = '../exps_withCori/'
% expname='topo4_H500Lx3k_s1.4dz1dx3n-20sm100_kv8e-5'
% expname = 'topo0_H500_s0.0016dz1dx3ln200n-20sm100_kv2e-4';
% expdir = '../exps_hires/';
% expname = 'hires_topo4_s0.0013_dz1dx6n-20';
% expname = 'topo4_H500_smo100m_s0.0014_dz1dx3ln200n-20'
% expdir = '../exps_topo4/';
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
% o = 12*26+3;
o = 12*20+3;

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
N2(:,1:Nr-1) = -cosd(topo)*gravity/rhoConst.*(rho(:,1:end-1)-rho(:,2:end))./(zz(1:end-1)-zz(2:end));

rho_ugrid = 0.5*(rho([Nx 1:Nx-1],:)+rho);
N2x = -sind(topo)*gravity/rhoConst.*(rho_ugrid([2:Nx 1],:)-rho_ugrid)./delX(1);
N2x_lower = NaN*zeros(Nx,Nr);
N2x_lower(:,1:Nr-1) = 0.5*(N2x(:,1:end-1)+N2x(:,2:end));
N2 = N2+N2x_lower;

% N2 = N2+N2const;

%--- make fig2


%%
figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 900 950]);

%%% coordinate
ax1 = subplot('position',[.03 .785 .3 .2]);
annotation('textbox',[0 0.993 0.15 0.01],'String','A','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
annotation('textbox',[0.07 0.993 0.3 0.01],'String','Slope-aligned coordinate','FontSize',fontsize+4,'interpreter','latex','LineStyle','None');
imshow('fig2/coordinate.png')

%--- Load MITgcm simulation
filename = [expdir expname '/RMSE_tt.mat'];
load(filename)
load('fig2/fig2_new_1.7.mat')
YLIM = [0 250];
% time_tidal=time_h/12;

%%% TKE time series
ax2 = subplot('position',[0.435 0.815 0.505 0.16]);
annotation('textbox',[0.38 0.993 0.15 0.01],'String','B','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
plot(time_tidal,log(pe)/2,'LineWidth',2);
hold on;
plot(time_tidal,log(ke)/2,'LineWidth',2);
plot(xxplot(fit_span)/12, y_fit(fit_span),'k--','LineWidth',2);
xlabel('Time (tidal cycles)','interpreter','latex');
ylabel('log(energy)','interpreter','latex');
h2 = legend('Turbulent potential energy','Turbulent kinetic energy','Linear fit',...
    'Fontsize',fontsize+1,'Position',[0.69 0.825 0.2285 0.0661],'interpreter','latex');
set(gca,'Fontsize',fontsize)
title('Normalized turbulent energy in the shear layer','Fontsize',fontsize+4,'interpreter','latex');
grid on;grid minor;
hold on;
ylim([-38 2])
xlim([0 28])


mycolor=cmocean('balance');
mycolor=mycolor(20:end-20,:);

%%% Velocity
ax3 = subplot('position',[0.07 0.62 0.87 0.125]);
annotation('textbox',[0 0.755 0.15 0.01],'String','C','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(time_tidal,zz-botZ,uu_timeseries');
hold on;shading interp;
contour(time_tidal,zz-botZ,uu_timeseries',[0.15:0.15:0.75],'color',darkgray)
contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1)
contour(time_tidal,zz-botZ,uu_timeseries',[-0.75:0.15:-0.15],'--','color',darkgray)
clim([-0.401 0.401])
ylabel('HAB (m)','interpreter','latex');
set(gca,'Fontsize',fontsize);
title('Across-isobath tidal velocity $u$','Fontsize',fontsize+4,'interpreter','latex','Position',[15,295-50])
ylim(YLIM)
h3=colorbar(ax3);
set(h3,'Position',[0.95    0.62    0.008    0.11]);
set(get(h3,'Title'),'String',{'$\ \ \ \ (\mathrm{m/s})$'},'interpreter','latex','FontSize',fontsize);
set(gca,'xtick',[])
colormap(mycolor);
freezeColors;
XLIM=[0 28];
xlim(XLIM);

%%% Temperature
ax4 = subplot('position',[0.07 0.46 0.87 0.13]);
annotation('textbox',[0 0.595 0.15 0.01],'String','D','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(time_tidal,zz-botZ,tt_timeseries');
hold on;shading interp;
contour(time_tidal,zz-botZ,uu_timeseries',[0.15:0.15:0.75],'color',darkgray)
contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1)
contour(time_tidal,zz-botZ,uu_timeseries',[-0.75:0.15:-0.15],'--','color',darkgray)
ylabel('HAB (m)','interpreter','latex');
clim([-0.1 0.1]);
set(gca,'Fontsize',fontsize);
title('Time-varying component of temperature','Fontsize',fontsize+4,'interpreter','latex','Position',[15,295-50])
ylim(YLIM)
h4=colorbar(ax4);
set(h4,'Position',[0.95    0.46   0.008    0.11]);
set(get(h4,'Title'),'String',{'$\ \ \ \ (^\circ \mathrm{C})$'},'interpreter','latex','FontSize',fontsize);
set(gca,'xtick',[])
colormap(mycolor);
freezeColors;
xlim(XLIM);

%%% Temperature snapshot
ax6 = subplot('position',[0.07 0.05 0.37 0.18]);
annotation('textbox',[0 0.24 0.15 0.01],'String','F','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(xx/1000,zz-botZ,tt');hold on;
shading interp;
% clim([-0.02 0.1])
clim([-0.02 0.06])
ylabel('HAB (m)','interpreter','latex');
set(gca,'Fontsize',fontsize);
ylim(YLIM)
h6=colorbar(ax6);
set(h6,'Position',[0.45    0.05   0.008    0.16]);
set(get(h6,'Title'),'String',{'$\ \ \ \ (^\circ \mathrm{C})$'},'interpreter','latex','FontSize',fontsize);
xlabel('$x$ (km)','interpreter','latex','FontSize',fontsize+2);
title('Temperature $T$ (snapshot)','Fontsize',fontsize+4,'interpreter','latex','Position',[0 295-50]);
% xlim([0 3])
colormap(mycolor);
freezeColors;

%%% N2 snapshot
ax7 = subplot('position',[0.57 0.05 0.37 0.18]);
annotation('textbox',[0.5 0.24 0.15 0.01],'String','G','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(xx/1000,zz-botZ,N2')
hold on;
contour(xx/1000,zz-botZ,N2',[0 0],'Color','c','LineWidth',1);
hold off;
shading interp;set(gca,'Fontsize',fontsize);
clim(([-1 1]+1)/1e6)
ylabel('HAB (m)','interpreter','latex');
set(gca,'Fontsize',fontsize);
ylim(YLIM)
h7=colorbar(ax7);
set(h7,'Position',[0.95    0.05   0.008    0.16]);
set(get(h7,'Title'),'String',{'$\ \ \ \ (1/\mathrm{s}^2)$',''},'interpreter','latex','FontSize',fontsize);
xlabel('$x$ (km)','interpreter','latex','FontSize',fontsize+2);
title('$\partial_{\tilde z} \mathcal B$ (snapshot)','Fontsize',fontsize+4,'interpreter','latex','Position',[0 297-50])
% xlim([0 3])
colormap(mycolor);
freezeColors;


%%% N2
ax5 = subplot('position',[0.07 0.3 0.87 0.13]);
annotation('textbox',[0 0.435 0.15 0.01],'String','E','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(time_tidal,zz-botZ,(N2_timeseries)')
hold on;
contour(time_tidal,zz-botZ,(N2_timeseries)',[0 0],'Color','c','LineWidth',1);
contour(time_tidal,zz-botZ,uu_timeseries',[0.15:0.15:0.75],'color',darkgray)
contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1)
contour(time_tidal,zz-botZ,uu_timeseries',[-0.75:0.15:-0.15],'--','color',darkgray)
shading interp;
xlabel('Time (tidal cycles)','interpreter','latex');
ylabel('HAB (m)','interpreter','latex');
set(gca,'Fontsize',fontsize);
title('Vertical buoyancy gradient $\partial_{\tilde z} \mathcal B$','Fontsize',fontsize+4,'interpreter','latex','Position',[15,297-50])
clim(([-1 1]+1)/1e6)
ylim(YLIM)
% colormap(cmocean('diff'));
h5=colorbar(ax5);
set(h5,'Position',[0.95    0.3   0.008    0.1]);
set(get(h5,'Title'),'String',{'$\ \ \ \ (1/\mathrm{s}^2)$',''},'interpreter','latex','FontSize',fontsize);
xlim(XLIM);

%%% Save the figure

print('-dpng','-r300','fig2/fig2.png');