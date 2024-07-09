clear;close all;
addpath ../analysis/colormaps/
fontsize = 15;
% load_colors;

addpath ../analysis/
addpath ../analysis/functions/
expname = 'topo0_H500_s0.0017dz1dx3ln200n-20sm100_kv2e-4';
expdir = '../exps_hires/';
load([expdir expname '/setParams.mat'])
loadexp;
rhoConst = 999.8;
fontsize = 16;

%-- Load initial temperature
fid = fopen(fullfile(exppath,'input','hydrogThetaFile.bin'),'r','b');
iniT = zeros(Nx,Ny,Nr);
for k=1:Nr
  iniT(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);
iniT = squeeze(iniT);


%-- Initial velocity


mycolor=cmocean('balance');
mycolor=mycolor(30:end-30,:);

%%
figure(1)
clf;
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 800 600]);

%--- Initial temperature field, indicating the sponge layer
ax1 = subplot('position',[0.1 0.78 0.8 0.17]);
annotation('textbox',[0 0.99 0.15 0.01],'String','a','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(xx,zz,iniT');shading flat;colormap(mycolor);
clim([-1.01 1.01]*1e-20)
ylim([-500 0])
yticks([-500 -250 0])
ylabel('$z$ (m)','interpreter','latex')
xlabel('$x$ (m)','interpreter','latex')
set(gca,'Fontsize',fontsize);
h1 = colorbar(ax1);
set(h1,'Position',[0.9050 0.7800 0.007 0.1700]);
set(get(h1,'Title'),'String','$(^\circ \mathrm{C})\ \ \ \ \ $','Fontsize',fontsize,'interpreter','latex');
title('Initial temperatue perturbation','Fontsize',fontsize+4,'interpreter','latex')

%--- Tidal amplitude and linear shear
ax2 = subplot('position',[0.07 0.05 0.4 0.58]);
annotation('textbox',[0 0.6 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');

set(gca,'Fontsize',fontsize);
title('Tidal amplitude and shear (linear)','Fontsize',fontsize+4,'interpreter','latex')

%--- Tidal amplitude and shear as a function of tanh(z)
ax3 = subplot('position',[0.57 0.05 0.4 0.58]);
annotation('textbox',[0.5 0.6 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');

set(gca,'Fontsize',fontsize);
title('Tidal amplitude and shear (tanh)','Fontsize',fontsize+4,'interpreter','latex')


% print('-dpng','-r200',['fig_supp/figS_gcm_configuration.png']);









% %--- Growth rate of the MITgcm simulation: linear shear, max tanh shear, mean tanh shear
% ax2 = subplot('position',[0.57 0.65 0.4 0.3]);
% annotation('textbox',[0.5 0.99 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
% set(gca,'Fontsize',fontsize);
% title('Growth rate (flat bottom)','Fontsize',fontsize+4,'interpreter','latex')
