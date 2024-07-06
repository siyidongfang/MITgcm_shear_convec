
% clear;close all;
addpath ../analysis/colormaps/
fontsize = 18;
load_colors;

%--- make fig3
load('fig3/MITgcm_growth_hires_topo4.mat')
load('fig3/topo4_kappa0.mat')

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 900 900]);

%%% coordinate
ax1 = subplot('position',[.048 .785 0.37 0.18]);
annotation('textbox',[0.028 0.993 0.15 0.01],'String','a','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
plot(shear_MITgcm,growth_MITgcm,'LineWidth',2);
grid on;axis tight;
hold on;
plot(shear_all,max_grow,'LineWidth',2);
ylabel('(hour$^{-1}$)','interpreter','latex');
set(gca,'Fontsize',fontsize);

