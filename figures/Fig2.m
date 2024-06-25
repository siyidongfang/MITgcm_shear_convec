
clear;close all;
addpath ../analysis/colormaps/
addpath freezeColors/
fontsize = 15;
load_colors;

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1000 950]);

%%% coordinate
ax1 = subplot('position',[.07 .795 .3 .2]);
imshow('fig2/coordinate.png')
annotation('textbox',[0 0.993 0.15 0.01],'String','A','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');

%--- Load MITgcm simulation
addpath ../analysis/
addpath ../analysis/functions/
expname = 'topo0_H500_s0.0016dz1dx3ln200n-20sm100_kv2e-4';
expdir = '../exps_hires/';
loadexp;
filename = [expdir expname '/RMSE.mat'];
load(filename)

%%% TKE time series
ax2 = subplot('position',[0.5 0.8 0.45 0.17]);
annotation('textbox',[0.45 0.993 0.15 0.01],'String','B','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');
plot(time_h/12,log(pe)/2,'LineWidth',2);
hold on;
plot(time_h/12,log(ke)/2,'LineWidth',2);
plot(xxplot(fit_span)/12, y_fit(fit_span),'k--','LineWidth',2);
xlabel('Time (tidal cycles)')
ylabel('log(energy)')
set(gca,'Fontsize',fontsize)
h2 = legend('Turbulent potential energy','Turbulent kinetic energy','Linear fit',...
    'Fontsize',fontsize+1,'Position',[0.7115 0.8146 0.2285 0.0661]);
title('Normalized turbulent energy in the shear layer','Fontsize',fontsize+3)
grid on;grid minor;
hold on;
ylim([-38 2])

%%% Velocity
ax3 = subplot('position',[0.1 0.6 0.8 0.17]);
annotation('textbox',[0.02 0.78 0.15 0.01],'String','C','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');

%%% Temperature

%%% N2

%%% Temperature snapshot

%%% N2 snapshot


%%% Save the figure

% print('-djpeg','-r300','fig2/fig2_matlab.png');
