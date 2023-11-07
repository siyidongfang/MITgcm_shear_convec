

clear;close all;
ne =1;
load_all

topo_slope = 4;
cos_slope = cosd(topo_slope);
sin_slope = sind(topo_slope);
N2_const = 1e-6;

No = nDumps;

%%% Predicted tidal velocity, include the velocity shear
utide = zeros(No,Nr);
% u_shear = 9e-4;
u_shear = 40e-4;
h_shear = 250;
amp_tide = u_shear * h_shear;
pi_const = 3.141592653589793;
omega_tide = 2*pi_const/43200;
shearProfile = ones(1,Nr);
for k=1:Nr
    if((zz(k)-zz(Nr))<h_shear)
        shearProfile(k) = (zz(k)-zz(Nr))/h_shear;
    end
end

for o=1:No
    nIter = dumpIters(o);
    time_s = nIter.*deltaT; %%% time in seconds
    utide(o,:) = amp_tide*cos(omega_tide*time_s)*shearProfile;
end

dz_uN_tide = zeros(1,Nr);
N2 = zeros(No,Nr);
N2(1,:) = N2_const;
dt = 3600;

for o=1:No
    uu_tide = squeeze(utide(o,:)); 
    uN_tide = uu_tide*N2_const*sin_slope;
    dz_uN_tide(2:Nr)=-diff(uN_tide)./delR(2:Nr);   %%% u-grid, upper level
    if(o>1)
        N2(o,:) = N2(o-1,:)-dz_uN_tide*dt;
    end
end

figure(1)
pcolor(N2');axis ij;shading flat;colorbar;colormap(redblue);clim([0 2e-6])


figure(3)
pcolor(utide');axis ij;shading flat;colorbar;colormap(redblue);clim([-0.2 0.2])

