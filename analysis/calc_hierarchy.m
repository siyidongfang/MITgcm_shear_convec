
clear;close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Calculate Critical Shear %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assume cross-canyon velocity is comparable to along-canyon velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topo_slope = 4;
topo_slope_radian = 4/180*pi;
cos_slope = cosd(topo_slope);
sin_slope = sind(topo_slope);
tan_slope = tand(topo_slope);
N2_const = 1e-6;
f0 = 1.18e-4;

pi_const = 3.141592653589793;
% Ttide = 43200;
% Ttide = 12.4206012*3600;
Ttide = 23.93447213*3600;
omega = 2*pi_const/Ttide;

Ric = 0.25; %%% Critical Richardson number
% Ric = 0.1;

Sc = omega/tan_slope;  %%% Critical shear predicted by the linear theory when N^2 receaches zero
d = N2_const*cos_slope/Ric; %%% Define a constant

Bu = N2_const*topo_slope_radian^2/f0^2; %%% Slope Burger Number

%%% Two solutions for the quadratic equation
S1 = -d/2/Sc - 1/2*sqrt((d/Sc)^2+4*d)
S2 = -d/2/Sc + 1/2*sqrt((d/Sc)^2+4*d)
S3 = d/2/Sc + 1/2*sqrt((d/Sc)^2+4*d)
S4 = d/2/Sc - 1/2*sqrt((d/Sc)^2+4*d)


4*omega^2*Ric*cos_slope/N2_const/sin_slope^2


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Calculate Critical Shear %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assume cross-canyon velocity is zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N2_const = 1e-6;
f0 = 1.18e-4;
pi_const = 3.141592653589793;
Ric = 0.25; %%% Critical Richardson number
tt = 1:60:86400*2;

topo_slope = 4;
Ttide = 23.93447213*3600;
omega = 2*pi_const/Ttide;
a = Ric*(cos(omega*tt)).^2;
b = N2_const*sind(topo_slope)/omega*sin(omega*tt);
c = -N2_const*cosd(topo_slope);

SS1_k1 = (-b+sqrt(b.^2-4*a*c))/2./a;
SS2_k1 = (-b-sqrt(b.^2-4*a*c))/2./a;

topo_slope = 4;
Ttide = 12.4206012*3600;
omega = 2*pi_const/Ttide;
a = Ric*(cos(omega*tt)).^2;
b = N2_const*sind(topo_slope)/omega*sin(omega*tt);
c = -N2_const*cosd(topo_slope);

SS1_m2 = (-b+sqrt(b.^2-4*a*c))/2./a;
SS2_m2 = (-b-sqrt(b.^2-4*a*c))/2./a;

topo_slope = 10;
Ttide = 12.4206012*3600;
omega = 2*pi_const/Ttide;
a = Ric*(cos(omega*tt)).^2;
b = N2_const*sind(topo_slope)/omega*sin(omega*tt);
c = -N2_const*cosd(topo_slope);

SS1_m2_slope6 = (-b+sqrt(b.^2-4*a*c))/2./a;
SS2_m2_slope6 = (-b-sqrt(b.^2-4*a*c))/2./a;

fontsize = 18;
figure(1)
clf;set(gcf,'color','w');
plot(tt/3600,SS1_k1,'LineWidth',3)
hold on;
plot(tt/3600,SS1_m2,'LineWidth',3)
plot(tt/3600,SS1_m2_slope6,'-.','LineWidth',3)
ylim([0.9 3]*1e-3);
hold on;
grid on;grid minor;
legend('Dinural tide, \theta_{topo}=4^o=0.07',...
    'Semidinural tide, \theta_{topo}=4^o=0.07',...
    'Semidinural tide, \theta_{topo}=6^o=0.10')
set(gca,'Fontsize',fontsize);
xlabel('Time (hours)')
ylabel('Critical Shear (1/s)')
title('The minimum shear required for Ri<0.25')





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Calculate Critical Shear %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assume cross-canyon velocity is 1/5 of the along-canyon velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N2_const = 1e-6;
f0 = 1.18e-4;
pi_const = 3.141592653589793;
Ric = 0.25; %%% Critical Richardson number
tt = 1:60:86400*2;

v_frac = 1/2; %%% cross-canyon velocity is 1/5 of the along-canyon velocity


topo_slope = 4;
Ttide = 23.93447213*3600;
omega = 2*pi_const/Ttide;
a = Ric*((cos(omega*tt)).^2+v_frac^2*(sin(omega*tt)).^2);
b = N2_const*sind(topo_slope)/omega*sin(omega*tt);
c = -N2_const*cosd(topo_slope);

SS1_k1 = (-b+sqrt(b.^2-4*a*c))/2./a;
SS2_k1 = (-b-sqrt(b.^2-4*a*c))/2./a;

topo_slope = 4;
Ttide = 12.4206012*3600;
omega = 2*pi_const/Ttide;
a = Ric*((cos(omega*tt)).^2+v_frac^2*(sin(omega*tt)).^2);
b = N2_const*sind(topo_slope)/omega*sin(omega*tt);
c = -N2_const*cosd(topo_slope);

SS1_m2 = (-b+sqrt(b.^2-4*a*c))/2./a;
SS2_m2 = (-b-sqrt(b.^2-4*a*c))/2./a;

topo_slope = 6;
Ttide = 12.4206012*3600;
omega = 2*pi_const/Ttide;
a = Ric*((cos(omega*tt)).^2+v_frac^2*(sin(omega*tt)).^2);
b = N2_const*sind(topo_slope)/omega*sin(omega*tt);
c = -N2_const*cosd(topo_slope);

SS1_m2_slope6 = (-b+sqrt(b.^2-4*a*c))/2./a;
SS2_m2_slope6 = (-b-sqrt(b.^2-4*a*c))/2./a;

fontsize = 18;
figure(1)
clf;set(gcf,'color','w');
plot(tt/3600,SS1_k1,'LineWidth',3)
hold on;
plot(tt/3600,SS1_m2,'LineWidth',3)
plot(tt/3600,SS1_m2_slope6,'-.','LineWidth',3)
ylim([0.9 3]*1e-3);
hold on;
grid on;grid minor;
legend('Dinural tide, \theta_{topo}=4^o=0.07',...
    'Semidinural tide, \theta_{topo}=4^o=0.07',...
    'Semidinural tide, \theta_{topo}=6^o=0.10')
set(gca,'Fontsize',fontsize);
xlabel('Time (hours)')
ylabel('Critical Shear (1/s)')
title('The minimum shear required for Ri<0.25')





%%

ne =1;
load_all

% No = nDumps;
No =24

%%% Predicted tidal velocity, include the velocity shear
utide = zeros(No,Nr);
% u_shear = 9e-4;
shearProfile = ones(1,Nr);
for k=1:Nr
    if((zz(k)-zz(Nr))<h_shear)
        shearProfile(k) = (zz(k)-zz(Nr))/h_shear;
    end
end

for o=1:No
    nIter = dumpIters(o);
    time_s = nIter.*deltaT; %%% time in seconds
    utide(o,:) = amp_tide*sin(omega*time_s+0.25*pi)*shearProfile;
end

dz_uN_tide = zeros(1,Nr);
N2 = zeros(No,Nr);
tt = zeros(No,Nr);

% N2(1,:) = N2_const;
dt = 3600;
hab = zz(end)-zz;

for o=1:No
    uu_tide = squeeze(utide(o,:)); 
    uN_tide = uu_tide*N2_const*sin_slope;
    dz_uN_tide(2:Nr)=-diff(uN_tide)./delR(2:Nr);   %%% u-grid, upper level
    if(o>1)
        N2(o,:) = N2(o-1,:)-dz_uN_tide*dt;
    end
    tt(o,:) = N2(o,:)/gravity/tAlpha.*hab;
end

figure(1)
pcolor(N2');axis ij;shading flat;colorbar;colormap(redblue);
% clim([0 2e-6])
clim([-1e-6 1e-6])

figure(2)
pcolor(tt');axis ij;shading flat;colorbar;colormap(redblue);
clim([-0.1 0.1])

figure(3)
pcolor(utide');axis ij;shading flat;colorbar;colormap(redblue);clim([-0.2 0.2])

