clear;close all;
addpath ../analysis/colormaps/
fontsize = 16;
load_colors;

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1400 700]);

%--- T/S diagram
ax1 = subplot('position',[0.048 0.58 0.25 0.38]);
annotation('textbox',[0.028 0.993 0.15 0.01],'String','A','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');

%--- salinity
ax2 = subplot('position',[0.38 0.58 0.25 0.38]);
annotation('textbox',[0.36 0.993 0.15 0.01],'String','B','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');

%--- N2
ax3 = subplot('position',[0.709 0.58 0.25 0.38]);
annotation('textbox',[0.69 0.993 0.15 0.01],'String','C','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');

%--- Time series of depth-averaged N2 and Ri of the large-scale flow
ax4 = subplot('position',[0.048 0.07 0.25 0.38]);
annotation('textbox',[0.028 0.482 0.15 0.01],'String','D','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');

%--- dbdz using observed u, without N2
ax5 = subplot('position',[0.38 0.07 0.25 0.38]);
annotation('textbox',[0.36 0.482 0.15 0.01],'String','E','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');

%--- dbdz using linear-fit u, without N2
ax6 = subplot('position',[0.709 0.07 0.25 0.38]);
annotation('textbox',[0.69 0.482 0.15 0.01],'String','F','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');


% print('-djpeg','-r300','fig_supp/figS1.png');
