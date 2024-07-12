clear;close all;
addpath ../analysis/colormaps/
fontsize = 15;
load_colors;

addpath ../analysis/
addpath ../analysis/functions/
expname = 'topo0_H500_s0.0017dz1dx3ln200n-20sm100_kv2e-4';
expdir = '../exps_hires/';
load([expdir expname '/setParams.mat'])
loadexp;
rhoConst = 999.8;
fontsize = 16;
load_colors;

%-- Load initial temperature
fid = fopen(fullfile(exppath,'input','hydrogThetaFile.bin'),'r','b');
iniT = zeros(Nx,Ny,Nr);
for k=1:Nr
  iniT(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);
iniT = squeeze(iniT);


h_shear = 250;
Shear = 1.4e-3;

%-- Initial velocity

% Linear
Hmax = 500;
for i=1:Nr
   if((zz(i)+Hmax)<h_shear) 
       shearProfile(i)=(zz(i)+Hmax)/h_shear;
   else
       shearProfile(i)=1.;
   end 
end 

u_linear = Shear*h_shear*shearProfile;
Nshear_smooth_half = 100;
Nsmooth_span = Nshear_smooth_half*2+1;
u_linear_smooth = smooth(u_linear,Nsmooth_span);

shear_linear_smooth = diff(u_linear_smooth)./diff(zz);
shear_linear = diff(u_linear)./diff(zz);
z_mid = 0.5*(zz(1:end-1)+zz(2:end));

% tanh, zero velocity at sea floor
Hmax = 800;
zz_tanh = -799.5:-0.5;
zz_tanh_mid = 0.5*(zz_tanh(1:end-1)+zz_tanh(2:end));

u_tanh_bot = h_shear*Shear *(1+ tanh( (zz_tanh+Hmax/2) / (h_shear/2) )) /2;
shear_tanh_bot = diff(u_tanh_bot)./diff(zz_tanh);

% tanh, zero velocity at center
Hmax = 800;
u_tanh_center = h_shear*Shear *(   tanh( (zz_tanh+Hmax/2) / (h_shear/2) )) /2;
shear_tanh_center = diff(u_tanh_center)./diff(zz_tanh);


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
xlim([-1500 1500])
ylabel('$z$ (m)','interpreter','latex')
xlabel('$x$ (m)','interpreter','latex')
set(gca,'Fontsize',fontsize);
h1 = colorbar(ax1);
set(h1,'Position',[0.9050 0.7800 0.007 0.1700]);
set(get(h1,'Title'),'String','$(^\circ \mathrm{C})\ \ \ \ \ $','Fontsize',fontsize,'interpreter','latex');
title('Initial temperatue perturbation','Fontsize',fontsize+4,'interpreter','latex')


%--- Tidal amplitude and linear shear
ax2 = axes('position',[0.08 0.08 0.39 0.53]);
annotation('textbox',[0.02 0.68 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
plot(ax2, u_linear_smooth, zz, 'k-', 'LineWidth', 2)
ylim(ax2, [-500 0])
hold(ax2, 'on');
plot(ax2, u_linear, zz, '--', 'LineWidth', 1.4, 'Color', green)
grid(ax2, 'on');
grid(ax2, 'minor');
hold(ax2, 'off');
set(ax2, 'FontSize', fontsize);
xlabel(ax2, 'Tidal amplitude, quasi-linear (m/s)', 'Fontsize', fontsize+4, 'interpreter', 'latex')
ylabel(ax2, '$z$ (m)', 'interpreter', 'latex')

% Secondary axes for shear
ax22 = axes('Position', get(ax2, 'Position'), 'Color', 'none');
plot(ax22, shear_linear_smooth/1e-3, z_mid, '-', 'Color', brown, 'LineWidth', 2);
set(ax22, 'YAxisLocation', 'right', 'XAxisLocation', 'top', 'Color', 'none', 'XColor', brown, 'FontSize', fontsize);
xlabel(ax22, 'Shear (10$^{-3}\,$s$^{-1}$)', 'Color', brown, 'Interpreter', 'latex', 'Fontsize', fontsize+4);
ylim(ax22, [-500 0])
xlim(ax22, [-0.25 1.75])
% xlim(ax22, [0 1.5])
ax22.YTick = [];

% Create a dummy axes to show bottom x-axis in black
ax23 = axes('Position', get(ax2, 'Position'), 'Color', 'none', 'XAxisLocation', 'bottom', 'YAxisLocation', 'right', 'Color', 'none', 'XColor', 'k', 'YColor', 'none');
ax23.XTick = ax2.XTick;  % Match the ticks
ax23.XTickLabel = [];    % Hide the labels
set(ax23, 'FontSize', fontsize);
xlim([0 0.4])

% Link the axes
linkaxes([ax2, ax22, ax23], 'y');



%--- Tidal amplitude and shear as a function of tanh(z)
ax3 = axes('position',[0.58 0.08 0.39 0.53]);
annotation('textbox',[0.52 0.68 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
plot(u_tanh_bot,zz_tanh,'k-','LineWidth',2)
hold on;
plot(u_tanh_center,zz_tanh,'k--','LineWidth',2)
set(gca,'Fontsize',fontsize);
grid on;grid minor;
xlim(ax3,[-0.2 0.4])
xlabel('Tidal amplitude, tanh (m/s)','Fontsize',fontsize+4,'interpreter','latex')
ylabel('$z$ (m)','interpreter','latex')


% Secondary axes for shear
ax32 = axes('Position', get(ax3, 'Position'), 'Color', 'none');
plot(ax32, shear_tanh_bot/1e-3, zz_tanh_mid, '-', 'Color', brown, 'LineWidth', 2);
% hold on;
% % plot(ax32, shear_tanh_center, zz_tanh_mid, '--', 'Color', brown, 'LineWidth', 2);
% hold off;
set(ax32, 'YAxisLocation', 'right', 'XAxisLocation', 'top', 'Color', 'none', 'XColor', brown, 'FontSize', fontsize);
xlabel(ax32, 'Shear (10$^{-3}\,$s$^{-1}$)', 'Color', brown, 'Interpreter', 'latex', 'Fontsize', fontsize+4);
ylim(ax32, [-800 0])
% xlim(ax32, [-0.5 2.5])
xlim(ax32, [-0.35 1.75])
 % xlim(ax32, [0 1.5])
% ax32.XTick = [-0.5:0.5:2] * 1e-3;
% ax32.XTickLabel = ['-0.5','0','-0.5','1','1.5','2'];
ax32.YTick = [];

% Create a dummy axes to show bottom x-axis in black
ax33 = axes('Position', get(ax3, 'Position'), 'Color', 'none', 'XAxisLocation', 'bottom', 'YAxisLocation', 'right', 'Color', 'none', 'XColor', 'k', 'YColor', 'none');
ax33.XTick = ax3.XTick;  % Match the ticks
ax33.XTickLabel = [];    % Hide the labels
set(ax33, 'FontSize', fontsize);
xlim(ax33,[-0.2 0.4])

% Link the axes
linkaxes([ax3, ax32, ax33], 'y');

print('-dpng','-r300',['fig_supp/figS_gcm_configuration_matlab.png']);





% %--- Growth rate of the MITgcm simulation: linear shear, max tanh shear, mean tanh shear
% ax2 = subplot('position',[0.57 0.65 0.4 0.3]);
% annotation('textbox',[0.5 0.99 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
% set(gca,'Fontsize',fontsize);
% title('Growth rate (flat bottom)','Fontsize',fontsize+4,'interpreter','latex')
