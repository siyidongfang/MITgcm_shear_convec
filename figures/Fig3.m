
clear;close all;
addpath ../analysis/colormaps/
fontsize = 18;
load_colors;

%--- make fig3

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 800 800]);

%%% coordinate
ax1 = subplot('position',[.03 .785 .3 .2]);
annotation('textbox',[0 0.993 0.15 0.01],'String','A','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');
