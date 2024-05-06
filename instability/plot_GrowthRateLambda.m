
clear; close all;
fontsize = 20;
load_colors;
addpath products/
load('GrowthRate-exps_test_Nr100.mat')
% load('GrowthRate_RK4_Nr200.mat')
% load('GrowthRate.mat')
% load('GrowthRate_new.mat')
% GrowthRate = squeeze(max(max(GrowthRate)))';
GrowthRate = max(GrowthRate_b2);
figure(1)
clf;
plot(Shear_parm,GrowthRate)

kx = 2*pi./lambda_parm;

figure(40)
clf;set(gcf,'Color','w')
zcolor=Shear_parm;
cmap = colormap(WhiteBlueGreenYellowRed(8));
zmap = linspace( min(zcolor), max(zcolor), length(cmap));
i=1;
color = interp1( zmap, cmap, zcolor(i));
semilogx(kx,GrowthRate_b2(:,i),'LineWidth',1.5,'Color',color);
grid on;grid minor;
hold on;
for i=2:length(Shear_parm)
    color = interp1( zmap, cmap, zcolor(i));
    semilogx(kx,GrowthRate_b2(:,i),'LineWidth',1.5,'Color',color);
end
hold off;
set(gca,'Fontsize',fontsize)
xlabel('Horizontal wavenumber, k_x (1/m)')
ylabel('Growth rate (1/hour)')
clim([min(zcolor), max(zcolor)])
c=colorbar;
grid on;

% load('GrowthRate_RK4_nu5e-6.mat')
% load('GrowthRate_lores_H150.mat')
% growthrate = squeeze(GrowthRate(4,:,:));
% growthrate = GrowthRate;
growthrate = GrowthRate_b2(:,1:10);
max_growth_allLambda = max(growthrate);

%%% Calculate the minimum Richardson Number

Shear_parm = [0.1:0.2:1.9]*1e-3;
Ri_min = NaN*zeros(1,length(Shear_parm));

NTtide = 10;
dt = 0.001;
Lt = NTtide*2*pi; % dimensionless simulation time
Nt = round(Lt/dt);
dt_real = NTtide*43200/Nt;
tt = dt:dt:Nt*dt;

Ptide = 43200;
omega = 2*pi/Ptide;
Nsquare = 1e-6;
topo = 4;    
for ns=1:length(Shear_parm)

    Shear = Shear_parm(ns);

    Ri_inverse = (Shear*cos(tt)).^2./(Nsquare*cosd(topo) - Nsquare*sind(topo)/omega*Shear*sin(tt));
    
    if(min(Ri_inverse)<0)
        isConvec(ns) = 1;
    end
    Ri_min(ns) = 1/max(Ri_inverse);

end


figure(7)
clf
set(gcf,'color','w','Position',[184 283 666 420]);
plot(Shear_parm,1./Ri_min,'LineWidth',2)
hold on;
plot(Shear_parm,1/0.25*ones(1,length(Shear_parm)),'--','LineWidth',2)
set(gca,'Fontsize',fontsize);
ylabel('Ri_{min}^{-1}')
xlabel('Shear \Lambda (1/s)')
grid on;




figure(11)
clf;
set(gcf,'color','w','Position',[184 283 666 420]);
plot(1./Ri_min,growthrate,'LineWidth',2);
hold on;
plot((1/0.25*ones(1,length(0:0.01:0.25))),0:0.01:0.25,'--','LineWidth',2)
grid on;grid minor;
set(gca,'Fontsize',fontsize);
xlabel('Ri_{min}^{-1}')
ylabel('Growth Rate (hour^{-1})')
% title('Maximum growth rate for all wavenumbers')
title('Growth rate (lambda = 100m)')



figure(12)
clf;
set(gcf,'color','w','Position',[184 283 666 420]);
plot(Shear_parm,growthrate,'LineWidth',2);
hold on;
% plot((1/0.25*ones(1,length(0:0.01:0.25))),0:0.01:0.25,'--','LineWidth',2)
grid on;grid minor;
set(gca,'Fontsize',fontsize);
xlabel('Shear \Lambda (1/s)')
ylabel('Growth Rate (hour^{-1})')
% title('Maximum growth rate for all wavenumbers')
title('Growth rate (lambda = 100m)')


%%

% figure(4)
% clf;
% set(gcf,'color','w','Position',[184 283 666 420]);
% pcolor(2*pi./(lambda_parm),log10(Ri_min),growthrate')
% shading flat;colorbar;
% colormap(WhiteBlueGreenYellowRed(0));
% set(gca,'Fontsize',fontsize);
% ylabel('log(Ri_{min})')
% xlabel('Along-slope wavenumber k (1/m)')
% title('Growth Rate (hour^{-1})')


figure(5)
clf;
set(gcf,'color','w','Position',[184 283 666 420]);
pcolor(log10(lambda_parm),1./(Ri_min),growthrate')
hold on;
plot(log10(lambda_parm),(1/0.25*ones(1,length(lambda_parm))),'--','Color',black,'LineWidth',3)
shading flat;colorbar;
colormap(WhiteBlueGreenYellowRed(0));
set(gca,'Fontsize',fontsize);
ylabel('Ri_{min}^{-1}')
xlabel('Along-slope wavelength log(\lambda) (m)')
title('Growth Rate (hour^{-1})')



figure(6)
clf;
set(gcf,'color','w','Position',[184 283 666 420]);
plot(1./(Ri_min),max_growth_allLambda,'LineWidth',2);
hold on;
plot(1/0.25*ones(1,length(0:0.01:0.55)),0:0.01:0.55,'--','LineWidth',2)
grid on;grid minor;
set(gca,'Fontsize',fontsize);
xlabel('Ri_{min}^{-1}')
ylabel('Growth Rate (hour^{-1})')
title('Maximum growth rate across all wavenumbers')





figure(8)
clf;
set(gcf,'color','w','Position',[184 283 666 420]);
plot(1./Ri_min(3:end),max_growth_allLambda(3:end),'LineWidth',2);
hold on;
plot((1/0.25*ones(1,length(0:0.01:0.55))),0:0.01:0.55,'--','LineWidth',2)
grid on;grid minor;
set(gca,'Fontsize',fontsize);
xlabel('Ri_{min}^{-1}')
ylabel('Growth Rate (hour^{-1})')
title('Maximum growth rate across all wavenumbers')




%%
% figure(1)
% clf;
% set(gcf,'color','w','Position',[133 306 717 582]);
% pcolor(2*pi./(lambda_parm),Shear_parm,growthrate')
% shading flat;colorbar;
% colormap(WhiteBlueGreenYellowRed(0));
% set(gca,'Fontsize',fontsize);
% ylabel('Shear \Lambda (1/s)')
% xlabel('Along-slope wavenumber k (1/m)')
% title('Growth Rate (hour^{-1})')


figure(2)
clf;
set(gcf,'color','w','Position',[184 283 666 420]);
pcolor(log10(lambda_parm),Shear_parm,growthrate')
hold on;
% plot(log10(lambda_parm),0.001*ones(1,length(lambda_parm)),'--','Color',red,'LineWidth',3)
shading flat;colorbar;
colormap(WhiteBlueGreenYellowRed(0));
set(gca,'Fontsize',fontsize);
ylabel('Shear \Lambda (1/s)')
xlabel('Along-slope wavelength log(\lambda) (m)')
title('Growth Rate (hour^{-1})')


figure(3)
clf;
set(gcf,'color','w','Position',[184 283 666 420]);
plot(Shear_parm,max_growth_allLambda,'LineWidth',2);
grid on;grid minor;
set(gca,'Fontsize',fontsize);
xlabel('Shear \Lambda (1/s)')
ylabel('Growth Rate (hour^{-1})')
title('Maximum growth rate across all wavenumbers')

