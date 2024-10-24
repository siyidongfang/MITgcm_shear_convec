clear all;close all

addpath ../observations/
addpath ../analysis/colormaps/

fontsize = 20;
load_colors;

load('MAVS2_velocity.mat')
zidx = 4:18; %%% Full depth
t_idx_plot = 1:1000;

CLIM = [-0.5 0.5];


fg1 = figure(1);
clf;
set(gcf,'Color','w','Position',[100 153 1100 600])
tiledlay = tiledlayout(2,1);

nexttile;
pcolor(time(t_idx_plot)/86400,depth(zidx),uu(zidx,t_idx_plot));
shading interp;colorbar;clim(CLIM)
ylabel('Depth (m)','Interpreter','latex');xlabel('Time (days)','Interpreter','latex');set(gca,'FontSize',fontsize)
axis ij;title('Along-canyon velocity at MAVS2 (m/s)','Interpreter','latex','FontSize',fontsize+5)


% colormap(redblue);
colormap(cmocean('balance'));


nexttile;
pcolor(time(t_idx_plot)/86400,depth(zidx),vv(zidx,t_idx_plot));
shading interp;colorbar;clim(CLIM)
ylabel('Depth (m)','Interpreter','latex');xlabel('Time (days)','Interpreter','latex');set(gca,'FontSize',fontsize)
axis ij;title('Cross-canyon velocity at MAVS2 (m/s)','Interpreter','latex','FontSize',fontsize+5)



tiledlay.TileSpacing = 'compact';
tiledlay.Padding = 'compact';

addpath ~/MITgcm_shear_convec/figures/
AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal')

print('-dpng','-r300','fig_supp_new/figS_obs_crosscanyon.png');
