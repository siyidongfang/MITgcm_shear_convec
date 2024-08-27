
clear;close all;

omega = 2*pi/43200;
N = 1e-3;

% inverseRi = 0:0.1:5;


load('../figures/fig4/Ri_flat.mat')
load('../figures/fig5/fig5_topo0_noDiff_eig.mat')
addpath ../analysis/colormaps/
fontsize = 18;
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_eig(i) = Ri_min(b(i));
end
max_grow_eig = max(grow_all);
rw_eig = lam_z_all./lam_x_all;


grow_smallm0 = grow_all(800,:);
grow_smallm0(grow_smallm0<0)=0;
grow_smallm0(isnan(grow_smallm0))=0;
figure(2)
plot(1./Ri_eig,grow_smallm0)

figure(1)
alpha = omega/N/2
pcolor(1./Ri_eig,1./rw_eig,grow_all);
hold on;
xx = 1./Ri_eig;
yy = (sqrt(alpha-omega^2/N^2)+sqrt(xx))*N/omega;
% yy = sqrt(xx+0.7)*N/omega;
% yy = sqrt(xx)*N/omega;
plot(xx,yy,'LineWidth',2)
% set(gca,'YScale','log')
grid on;
grid minor;
hold on;
shading interp;
set(gca,'Fontsize',fontsize);
colormap(WhiteBlueGreenYellowRed(0))
colorbar;
clim([0 0.4])
xlim([0 5])