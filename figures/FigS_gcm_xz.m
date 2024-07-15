
clear;
close all;
addpath ../analysis/colormaps/
fontsize = 15;
load_colors;

addpath ../analysis/
addpath ../analysis/functions/
expname = 'topo4_H500_smo100m_s0.0014_dz1dx3ln200n-20'
expdir = '../exps_topo4/';
loadexp;
rhoConst = 999.8;
YLIM = [0 400];

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(1)); 
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

botZ =zz(end);

%--- snapshots
o = 12*26+3;
tt = squeeze(rdmds([exppath,'/results/THETA'],dumpIters(o)));
uu = squeeze(rdmds([exppath,'/results/UVEL'],dumpIters(o)));
vv = squeeze(rdmds([exppath,'/results/VVEL'],dumpIters(o)));
ww = squeeze(rdmds([exppath,'/results/WVEL'],dumpIters(o)));

tt(tt==0)=NaN;
ww(ww==0)=NaN;
uu(uu==0)=NaN;
vv(vv==0)=NaN;

S = NaN*zeros(Nx,Nr);
S(:,1:Nr-1) = (uu(:,1:end-1)-uu(:,2:end))./(zz(1:end-1)-zz(2:end));

rho = rhoConst.*(1-(tt-tRef)*tAlpha);
N2 = NaN*zeros(Nx,Nr);
N2(:,1:Nr-1) = -gravity/rhoConst.*(rho(:,1:end-1)-rho(:,2:end))./(zz(1:end-1)-zz(2:end));

Ri = NaN*zeros(Nx,Nr);
Ri = N2./(S.^2);

Hz = sum(delR);
N2const = (1e-3)^2;
tNorth = N2const *(zz+Hz) /9.81/2e-4;
tt_background = ones(Nx,Nr);
for k=1:Nr
    tt_background(:,k) = squeeze(tt_background(:,k))*tNorth(k);
end
tt = tt + tt_background;


dz = delR(1);
dx = delX(1);
du_dz = zeros(Nx,Nr);
dv_dz = zeros(Nx,Nr);
dw_dz = zeros(Nx,Nr);
du_dx = zeros(Nx,Nr);
dv_dx = zeros(Nx,Nr);
dw_dx = zeros(Nx,Nr);
w_tlevel = zeros(Nx,Nr);
du_dx_ugrid_wlev = zeros(Nx,Nr);
dv_dx_wlev = zeros(Nx,Nr);

du_dz(:,2:Nr) = - diff(uu,1,2) ./dz; 
dv_dz(:,2:Nr) = - diff(vv,1,2) ./dz; % on v-grid
dv_dz_ugrid = 0.5.*(dv_dz([Nx 1:Nx-1],:,:)+dv_dz(1:Nx,:,:));
w_tlevel(:,1:Nr-1) = 0.5*(ww(:,1:Nr-1)+ww(:,2:Nr));
dw_dz(:,2:Nr) = - diff(w_tlevel,1,2) ./dz; % on v-grid
dw_dz_ugrid = 0.5.*(dw_dz([Nx 1:Nx-1],:,:)+dw_dz(1:Nx,:,:));

du_dx= (uu([2:Nx 1],:)-uu) ./dx; % on v-grid
du_dx_ugrid = 0.5.*(du_dx([Nx 1:Nx-1],:,:)+du_dx(1:Nx,:,:));
du_dx_ugrid_wlev(:,2:Nr) = 0.5*(du_dx_ugrid(:,1:Nr-1)+du_dx_ugrid(:,2:Nr));

dw_dx= (ww-ww([Nx 1:Nx-1],:)) ./dx; 

dv_dx= (vv-vv([Nx 1:Nx-1],:)) ./dx; 
dv_dx_wlev(:,2:Nr)= 0.5*(dv_dx(:,1:Nr-1)+dv_dx(:,2:Nr));


epsilon = viscAr.*(du_dz.^2+dv_dz_ugrid.^2+dw_dz_ugrid.^2)  ...
    +viscAh.*(du_dx_ugrid_wlev.^2+dv_dx_wlev.^2+dw_dx.^2);  %%% Mass-point, lower level
    


% drho_dx = zeros(Nx,Nr);
% drho_dx(2:Nx,:)=diff(rho,1,1)./delX(1);
% db_dx = -gravity./rhoConst.*drho_dx;
% 
% drho_dz = zeros(Nx,Nr);
% drho_dz(:,1:Nr-1) = - diff(rho,1,2) ./ dz;
% db_dz = -gravity./rhoConst.*drho_dz;
% 
% chi_x = diffKhT.*(db_dx.^2); %%% Mass-point, cell center
% chi_z = diffKrT.*db_dz.^2;   %%% Mass-point, lower level
% chi_x_lower = zeros(Nx,Nr);
% chi_x_lower(:,1:Nr-1) = 0.5 * (chi_x(:,1:Nr-1) + chi_x(:,2:Nr));
% 
% chi_buoy = chi_z + chi_x_lower; %%% Mass-point, lower level


dT_dx= (tt-tt([Nx 1:Nx-1],:)) ./dx; % u-grid, t-lev

tt_ugrid_wlev = zeros(Nx,Nr);
tt_ugrid = 0.5*(tt([Nx 1:Nx-1],:)+tt);
tt_ugrid_wlev(:,2:Nr) =  0.5*(tt_ugrid(:,1:Nr-1)+tt_ugrid(:,2:Nr));

dT_dz = zeros(Nx,Nr);
dT_dz(:,1:Nr-1) = - diff(tt_ugrid_wlev,1,2) ./ dz;

chi_x = diffKhT.*(dT_dx.^2); %%% u-grid, t-lev
chi_z = diffKrT.*dT_dz.^2;   %%% u-grid, t-lev
chi_tt = chi_z + chi_x; %%% u-grid, t-lev


figure(1)
clf;
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 800 600]);



ax5 = subplot('position',[0.07 0.08 0.365 0.22]);
annotation('textbox',[0 0.345 0.15 0.01],'String','e','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(xx/1000,zz-botZ,log10(epsilon)');
shading flat;
ylabel('HAB (m)','interpreter','latex');
xlabel('$x$ (km)','interpreter','latex');
set(gca,'Fontsize',fontsize);
clim([-11 -8])
% clim([-11 -9.4])
% clim([0 1]*1e-9)
% title('Dissipation rate of TKE','Fontsize',fontsize+4,'interpreter','latex');
title('$\epsilon=\nu(\nabla \mathbf{u})^2$','Fontsize',fontsize+4,'interpreter','latex');
h5 = colorbar(ax5,'Position',[0.44 0.08+0.015 0.008 0.18]);
set(get(h5,'Title'),'String',{'$\log(\mathrm{m^2/s^3})$'},'interpreter','latex','FontSize',fontsize,'Position',[3.2000 128 0]);
colormap(cmocean('haline'))
freezeColors(ax5);
ylim(YLIM)


ax6 = subplot('position',[0.57 0.08 0.365 0.22]);
annotation('textbox',[0.5 0.345 0.15 0.01],'String','f','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(xx/1000,zz-botZ,log10(chi_tt)');
shading flat;
ylabel('HAB (m)','interpreter','latex');
xlabel('$x$ (km)','interpreter','latex');
set(gca,'Fontsize',fontsize);
clim([-12.5 -8.6])
% clim([-11 -8])
% clim([0 1]*1e-9)
% title('Dissipation rate of temp. variance','Fontsize',fontsize+4,'interpreter','latex');
title('$\chi=\kappa(\nabla T)^2$','Fontsize',fontsize+4,'interpreter','latex');
h6 = colorbar(ax6,'Position',[0.94 0.08+0.015 0.008 0.18]);
set(get(h6,'Title'),'String',{'$\log(\mathrm{^\circ C^2/s})$'},'interpreter','latex','FontSize',fontsize,'Position',[3.2000 128 0]);
colormap(cmocean('thermal'))
freezeColors(ax6);
ylim(YLIM)




ax1 = subplot('position',[0.07 0.72 0.365 0.22]);
annotation('textbox',[0 0.99 0.15 0.01],'String','a','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(xx/1000,zz-botZ,uu');
shading flat;
clim([-1 1]*0.18)
ylabel('HAB (m)','interpreter','latex');
set(gca,'Fontsize',fontsize);
title('Across-isobath velocity $u$','Fontsize',fontsize+4,'interpreter','latex');
colormap(cmocean('balance'))
h1 = colorbar(ax1,'Position',[0.44 0.72+0.015 0.008 0.18]);
set(get(h1,'Title'),'String',{'$\ \ \ \ (\mathrm{m/s})$'},'interpreter','latex','FontSize',fontsize,'Position',[3.2000 128 0]);
xlim([-1.5 1.5])
freezeColors(ax1);
ylim(YLIM)


ax2 = subplot('position',[0.57 0.72 0.365 0.22]);
annotation('textbox',[0.5 0.99 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(xx/1000,zz-botZ,tt');
shading flat;
clim([-0.01 0.21])
ylabel('HAB (m)','interpreter','latex');
set(gca,'Fontsize',fontsize);
title('Temperature $T$','Fontsize',fontsize+4,'interpreter','latex');
h2 = colorbar(ax2,'Position',[0.94 0.7200+0.015 0.008 0.18]);
set(get(h2,'Title'),'String',{'$\ \ \ \ (^\circ\mathrm{C})$'},'interpreter','latex','FontSize',fontsize,'Position',[3.2000 128 0]);
xlim([-1.5 1.5])
colormap(cmocean('balance'))
freezeColors(ax2);
ylim(YLIM)


ax3 = subplot('position',[0.07 0.4 0.365 0.22]);
annotation('textbox',[0 0.665 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(xx/1000,zz-botZ,S');
shading flat;
clim([0 0.7]*3e-3)
ylabel('HAB (m)','interpreter','latex');
set(gca,'Fontsize',fontsize);
title('Absolute shear $\vert\Lambda\vert=\vert\partial_z u\vert$','Fontsize',fontsize+4,'interpreter','latex');
h3 = colorbar(ax3,'Position',[0.44 0.4+0.015 0.008 0.18]);
set(get(h3,'Title'),'String',{'$\ \ \ \ (\mathrm{1/s})$'},'interpreter','latex','FontSize',fontsize,'Position',[3.2000 128 0]);
colormap(cmocean('haline'))
freezeColors(ax3);
ylim(YLIM)


ax4 = subplot('position',[0.57 0.4 0.365 0.22]);
annotation('textbox',[0.5 0.665 0.15 0.01],'String','d','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(xx/1000,zz-botZ,ww');
shading flat;
clim([-1 1]*0.03)
ylabel('HAB (m)','interpreter','latex');
set(gca,'Fontsize',fontsize);
title('Slope-normal velocity $w$','Fontsize',fontsize+4,'interpreter','latex');
h4 = colorbar(ax4,'Position',[0.94 0.4+0.015 0.008 0.18]);
set(get(h4,'Title'),'String',{'$\ \ \ \ (\mathrm{m/s})$'},'interpreter','latex','FontSize',fontsize,'Position',[3.2000 128 0]);
xlim([-1.5 1.5])
colormap(cmocean('balance'))
freezeColors(ax4);
ylim(YLIM)



print('-dpng','-r300',['fig_supp/figS_gcm_xz_matlab.png']);



figure(2)
clf;
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 800 600]);
ax3 = subplot('position',[0.07 0.4 0.365 0.22]);
clim([0 0.7]*3e-3)
set(gca,'Fontsize',fontsize);
h3 = colorbar(ax3,'Position',[0.44 0.4+0.015 0.008 0.18]);
ax5 = subplot('position',[0.07 0.08 0.365 0.22]);
clim([-11 -8])
% clim([-11 -9.4])
h5 = colorbar(ax5,'Position',[0.44 0.08+0.015 0.008 0.18]);
colormap(cmocean('haline'))
set(gca,'Fontsize',fontsize);
print('-dpng','-r300',['fig_supp/figS_gcm_xz_colorbar1.png']);



figure(3)
clf;
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 800 600]);
ax6 = subplot('position',[0.57 0.08 0.365 0.22]);
set(gca,'Fontsize',fontsize);
clim([-12.5 -8.6])
% clim([-11 -8])
h6 = colorbar(ax6,'Position',[0.94 0.08+0.015 0.008 0.18]);
colormap(cmocean('thermal'))
print('-dpng','-r300',['fig_supp/figS_gcm_xz_colorbar2.png']);




