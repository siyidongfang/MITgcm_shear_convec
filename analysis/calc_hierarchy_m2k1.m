
clear;close all;
fontsize = 18;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Calculate Critical Shear %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assume cross-canyon velocity is 1/5 of the along-canyon velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N2_const = 1e-6;
f0 = 1.18e-4;
pi = 3.141592653589793;
Ric = 0.25; %%% Critical Richardson number
tt = 1:60*30:86400*180;

v_frac = 0; %%% cross-canyon velocity is 1/5 of the along-canyon velocity

topo = 4;
T2 = 23.93447213*3600;
T1 = 12.4206012*3600;

om1 = 2*pi/T1;
om2 = 2*pi/T2;
phi2 = [0:1/1000:1]*2*pi;
r = [0:1/100:1];
min_S1c = NaN*zeros(length(r),length(phi2));
for k = 1:length(r)
    k
    for o = 1:length(phi2)
        clear S1c_tt
        S1c_tt = 1/tand(topo)./(1/om1*sin(om1*tt)+r(k)/om2*sin(om2*tt+phi2(o)));
        min_S1c(k,o) = min(abs(S1c_tt));
    end
end

% rr = 0.5
rr = 0.5
% pphi = 0.1*pi
pphi = 0
r0_S1c = 1/tand(topo)./(1/om1*sin(om1*tt)+rr/om2*sin(om2*tt+pphi));
figure(4)
clf;set(gcf,'color','w');
plot(tt(1:500)/86400,abs(r0_S1c(1:500)),'LineWidth',2);grid on;grid minor;
xlabel('Time, t (days)')
ylim([0 5]/1e3)
xlim([0 max(tt(1:500)/86400)])
ylabel('Critical Shear (1/s)')
set(gca,'Fontsize',fontsize);
title('Minimum shear required for $N^2<0$ (r=0.5, $\varphi=0$)','Interpreter','latex')
% title('Minimum shear required for $N^2<0$ (r=0)','Interpreter','latex')




S1c = mean(min_S1c,2);

for k = 1:length(r)
    E(k) = rmse(min_S1c(k,:),S1c(k)*ones(1,length(phi2)));
end

figure(1)
clf;set(gcf,'color','w');
plot(r,S1c,'LineWidth',2);
grid on;grid minor;
title('Critical shear for N^2<0')
ylabel('Critical Shear (1/s)')
xlabel('Ratio between K1 and M2 tidal current amplitudes, r')
set(gca,'Fontsize',fontsize);


figure(2)
clf;set(gcf,'color','w');
plot(r,E);
grid on;grid minor;
title('Root-mean-square error of the Critical Shear')
ylabel('RMSE')
xlabel('Ratio between K1 and M2 tidal current amplitudes, r')
set(gca,'Fontsize',fontsize);


figure(3)
clf;set(gcf,'color','w');
pcolor(r,phi2/pi,min_S1c');shading flat;colorbar;
title('Critical shear for N^2<0')
ylabel('$\varphi\ (\pi)$','Interpreter','latex')
xlabel('Ratio between K1 and M2 tidal current amplitudes, r')
set(gca,'Fontsize',fontsize);





% plot(tt(1:1000),abs(S1c(1:1000)));ylim([0 0.003])
% plot(phi2,min_S1c);grid on;grid minor;

%%

a = Ric*((cos(om2*tt)).^2+v_frac^2*(sin(om2*tt)).^2);
b = N2_const*sind(topo)/om2*sin(om2*tt);
c = -N2_const*cosd(topo);

SS1_k1 = (-b+sqrt(b.^2-4*a*c))/2./a;
SS2_k1 = (-b-sqrt(b.^2-4*a*c))/2./a;

topo = 4;
T2 = 12.4206012*3600;
om2 = 2*pi/T2;
a = Ric*((cos(om2*tt)).^2+v_frac^2*(sin(om2*tt)).^2);
b = N2_const*sind(topo)/om2*sin(om2*tt);
c = -N2_const*cosd(topo);

SS1_m2 = (-b+sqrt(b.^2-4*a*c))/2./a;
SS2_m2 = (-b-sqrt(b.^2-4*a*c))/2./a;

topo = 6;
T2 = 12.4206012*3600;
om2 = 2*pi/T2;
a = Ric*((cos(om2*tt)).^2+v_frac^2*(sin(om2*tt)).^2);
b = N2_const*sind(topo)/om2*sin(om2*tt);
c = -N2_const*cosd(topo);

SS1_m2_slope6 = (-b+sqrt(b.^2-4*a*c))/2./a;
SS2_m2_slope6 = (-b-sqrt(b.^2-4*a*c))/2./a;



% fontsize = 18;
% figure(1)
% clf;set(gcf,'color','w');
% plot(tt/3600,SS1_k1,'LineWidth',3)
% hold on;
% plot(tt/3600,SS1_m2,'LineWidth',3)
% plot(tt/3600,SS1_m2_slope6,'-.','LineWidth',3)
% ylim([0.9 3]*1e-3);
% hold on;
% grid on;grid minor;
% legend('Dinural tide, \theta_{topo}=4^o=0.07',...
%     'Semidinural tide, \theta_{topo}=4^o=0.07',...
%     'Semidinural tide, \theta_{topo}=6^o=0.10')
% set(gca,'Fontsize',fontsize);
% xlabel('Time (hours)')
% ylabel('Critical Shear (1/s)')
% title('The minimum shear required for Ri<0.25')





