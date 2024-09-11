
clear;close all;
addpath ../analysis/colormaps/
fontsize = 25;
load_colors;

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 950 900]);


load('fig4/fig4_km.mat')
load('../instability_km/budget/flat_output.mat')

% load('../instability_km/budget/topo4_output.mat')


tidx = nt_percycle*5+1:nt_percycle*7;
%- Calculate the buoyancy budget
uB0x = -re_uuu(tidx)*N^2*ss;
wB0z = -re_www(tidx)*N^2*cs;
wBz  =  re_www(tidx)*shear/omega*N^2*ss.*st(tidx);
% diffusion = 0
dbdt = [0 (re_buoy(3:end)-re_buoy(1:end-2))/dt/2 0];
dbdt = dbdt(tidx);

%- Normalization
uB0x = uB0x/max(abs(dbdt));
wB0z = wB0z/max(abs(dbdt));
wBz = wBz/max(abs(dbdt));
dbdt = dbdt/max(abs(dbdt));
tt = tt(tidx);


%--- Compute the vorticity budget
re_zeta = real(zeta);
dzetadt = [0 (re_zeta(3:end)-re_zeta(1:end-2))/dt/2 0];
dzetadt = dzetadt(tidx);

bxcs = real(1i*kx.*buoy(tidx)*cs);
bzss = real(-1i*mz_t(tidx).*buoy(tidx)*ss);

%--- Normalization
bxcs = bxcs/max(abs(dzetadt));
bzss = bzss/max(abs(dzetadt));
dzetadt = dzetadt/max(abs(dzetadt));


ax5 = subplot('position',[.06 .06+0.1 0.4 0.225*3]);
% annotation('textbox',[0.028+0.02 0.31 0.15 0.01],'String','e','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
lzres = plot(tt/43200,dzetadt-bxcs-bzss,'-','LineWidth',3,'Color',boxcolor);
hold on;
lbxcs = plot(tt/43200,bxcs,'LineWidth',2,'Color',blue);
lbzss = plot(tt/43200,bzss,'LineWidth',2,'Color',yellow);
ldzetadt = plot(tt/43200,dzetadt,'LineWidth',1,'Color',green);
grid on;grid minor;
xlabel('Time (tidal cycles)','interpreter','latex')
set(gca,'Fontsize',fontsize);
% title('Vorticity budget','interpreter','latex','Fontsize',fontsize+5);
title('Normalized vorticity budget','interpreter','latex','Fontsize',fontsize+5);
% ylim([-1 1])
l51 = legend([lbxcs,lbzss], ...
    '$\partial_x b^\prime \cos\theta$','-$\partial_z b^\prime \sin\theta$',...
    'interpreter','latex','Position',[0.075 0.225 0.0889 0.0668/4*3]);
ah=axes('position',get(ax5,'position'),'visible','off');
set(gca,'Fontsize',fontsize);
l52 = legend(ah,[ldzetadt lzres], ...
    '$\partial_\tau \zeta = \partial_t \zeta+ U\partial_x \zeta$','Residual',...
    'interpreter','latex','Position', [0.065 0.0825 0.1626 0.0458]');
legend('boxoff') 



ax6 = subplot('position',[.57 .06+0.1 0.4 0.225*3]);
% annotation('textbox',[0.538+0.02 0.31 0.15 0.01],'String','f','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
lres = plot(tt/43200,dbdt-uB0x-wB0z-wBz,'-','LineWidth',3,'Color',boxcolor);
hold on;
luB0x = plot(tt/43200,uB0x,'LineWidth',2,'Color',RED1);
lwBz = plot(tt/43200,wBz,'LineWidth',2,'Color',yellow);
lwB0z = plot(tt/43200,wB0z,'LineWidth',2,'Color',blue);
ldbdt = plot(tt/43200,dbdt,'LineWidth',1,'Color',green);
grid on;grid minor;
xlabel('Time (tidal cycles)','interpreter','latex')
set(gca,'Fontsize',fontsize);
% title('Buoyancy budget','interpreter','latex','Fontsize',fontsize+5);
title('Normalized buoyancy budget','interpreter','latex','Fontsize',fontsize+5);
% ylim([-1 1])
l61 = legend([lwB0z,lwBz,luB0x], ...
    '$w^\prime\partial_z B_0$','$w^\prime\partial_z B$','$u^\prime\partial_x B_0$',...
    'interpreter','latex','Position',[0.58 0.2104 0.0889 0.0668]);
ah=axes('position',get(ax6,'position'),'visible','off');
set(gca,'Fontsize',fontsize);
l62 = legend(ah,[ldbdt lres], ...
    '$\partial_\tau b^\prime = \partial_t b^\prime + U\partial_x b^\prime$','Residual',...
    'interpreter','latex','Position', [0.58 0.0825 0.1626 0.0458]');
legend('boxoff') 


figure(2)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 950 900]);

%--- plot vertical velocity and buoyancy perturbation
ax5 = subplot('position',[.06 .06+0.1 0.4 0.225*3]);
% annotation('textbox',[0.028+0.02 0.31 0.15 0.01],'String','e','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
plot(tt/43200,re_buoy(tidx)./max(abs(re_buoy(tidx))),'LineWidth',2,'Color',blue);
hold on;
plot(tt/43200,re_www(tidx)./max(abs(re_www(tidx))),'LineWidth',2,'Color',black);
plot(tt/43200,re_zeta(tidx)./max(abs(re_zeta(tidx))),'LineWidth',2,'Color',green);
grid on;grid minor;
l5 = legend('Buoyancy perturbation','Vertical velocity','Horizontal vorticity $\zeta$','interpreter','latex','Position',[0.07 0.2265 0.2063 0.0458]);
% l5 = legend('Buoyancy perturbation','Vertical velocity','interpreter','latex','Position',[0.07 0.2265 0.2063 0.0458]);
xlabel('Time (tidal cycles)','interpreter','latex')
set(gca,'Fontsize',fontsize);
title('Normalized perturbations','interpreter','latex','Fontsize',fontsize+5);
% title('Normalized perturbations','interpreter','latex','Fontsize',fontsize+5);

