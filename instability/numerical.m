addpath ../analysis/colormaps/

%%%% Define constants
Hshear = 300;
S = 0.8e-3; 
U0 = Hshear * S;
N = 1e-3;
topo = 4;
omega = 2*pi/43200;
nu = 2e-6;
kappa = 2e-6;
Pr = nu/kappa;
delta = sqrt(2*nu/omega);
Re = U0*delta/nu;
C = N*sind(topo)/omega;
D = Re*delta/2/Hshear;

%%% Model dimension
Lz = Hshear/delta;  % dimensionless domain height
dz = 10;
Nz = round(Lz/dz);
zz = dz:dz:Nz*dz;       % Height above topography
zzd = zz*delta;    % dimensional veriable

Lt = 1*86400*omega; % dimensionless simulation time
dt = 0.001;
Nt = round(Lt/dt);
tt = dt:dt:Nt*dt;
ttd = tt/omega;      % dimensional veriable

%%%% Define variables
psi = zeros(Nt,Nz);
zeta = zeros(Nt,Nz);
buoy = zeros(Nt,Nz);

dbdt = zeros(1,Nz);
dzetadt = zeros(1,Nz);

dbdz = zeros(1,Nz);
d2bdz2 = zeros(1,Nz);
dpsidz = zeros(1,Nz);
d2psidz2 = zeros(1,Nz);

kx = 0.01; %%% wave number

%%% Initial condition
psi(1,:) = 0;
zeta(1,:) = 0;
buoy(1,:) = 0;

%%%% Integrate non-dimensionalized b and zeta with time
%%%% 1st-order and 2nd-order centered difference
%%%% Euler forward scheme for time advancement
% %%%%  TO DO: Fourth-order Runge-Kutta scheme for time advancement

%%%% Ignore diffusion

for o=1:Nt-1

    psi(o,:) = cumsum(cumsum(zeta(o,:)*dz)*dz);

    for m = 2:Nz-1
        dpsidz(m)   = (psi(o,m+1)-psi(o,m-1))/dz;
        d2psidz2(m) = (psi(o,m-1)-2*psi(o,m)+psi(o,m+1))/dz^2;
        dbdz(m)   = (buoy(o,m+1)-buoy(o,m-1))/dz;
        d2bdz2(m) = (buoy(o,m-1)-2*buoy(o,m)+buoy(o,m+1))/dz^2;
    end

    dzetadt = 1i*kx*D*(dpsidz-zz.*(d2psidz2-kx^2*psi(o,:)))*cos(o*dt) ...
         +C^2*(1i*kx*cotd(topo)-dbdz);

    dbdt = dpsidz + (-1i*kx*cotd(topo)+1i*kx*D*sin(o*dt))*psi(o,:) ...
         -1i*kx*D*zz*cos(o*dt).*buoy(o,:);
         + (d2bdz2 - kx^2.*buoy(o,:))/2/Pr; % diffusion

    %%% Euler forward
    buoy(o+1,:) = dbdt*dt;
    zeta(o+1,:) = dzetadt*dt;

end

re_psi = real(psi); re_psi(re_psi==0)=NaN;
re_zeta = real(zeta); re_zeta(re_zeta==0)=NaN;
re_buoy = real(buoy); re_buoy(re_buoy==0)=NaN;


%%%%% Floquet stability




%%

figure(1)
fontsize = 16;
clf;set(gcf,'color','w');
subplot(3,1,1)
pcolor(ttd/3600,zzd,re_psi');shading flat;colorbar;colormap(redblue);
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Streamfunction \psi','Fontsize',fontsize+3);
clim([-100 100])

subplot(3,1,2)
pcolor(ttd/3600,zzd,re_zeta');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Horizontal vorticity perturbation \zeta','Fontsize',fontsize+3);
clim([-4 4]/1e4)


subplot(3,1,3)
pcolor(ttd/3600,zzd,re_buoy');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Buoyancy perturbation b','Fontsize',fontsize+3);
clim([0 0.01])








