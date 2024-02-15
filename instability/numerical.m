

outputname = [expdir 'output.mat'];

fontsize = 16;

%%%% Define constants
if(NOdiffusion)
    nu = 0;
    kappa = 0;
else
    nu = 2e-4; %%% Kaiser and Pratt 2022: nu=kappa=2e-6
    kappa = 2e-4;
end
Pr = nu/kappa;
t1hour = 3600;
omega = 2*pi/Ptide;

%%% Model dimension
Lz = Hdepth/delta;     % dimensionless domain height
Nr = round(Lz/dz);
zz = dz/2:dz:(Nr*dz-dz/2);  % Height above topography   
zz_wgrid = 0:dz:((Nr)*dz);
dz_real  = dz*delta;

Hshear = 300;
Lshear = Hshear/delta; % dimensionless vertical scale for velocity shear
Nshear = round(Lshear/dz);
% Hshear = zz(Nshear)*delta;
U0 = Uconst + Hshear * Shear;
Re = U0*delta/nu;

Lt = NTtide*43200*omega; % dimensionless simulation time
Nt = round(Lt/dt)
dt_real = NTtide*43200/Nt;
tt = dt:dt:Nt*dt;

C1 = U0/delta/omega;
C2 = N*sind(topo)/omega;
C3 = nu/delta^2/omega;
C4 = kappa/delta^2/omega;

%%%% Define variables
psi = zeros(Nt,Nr+1);
zeta = zeros(Nt,Nr+1);
buoy = zeros(Nt,Nr);

p0 = zeros(1,Nr+1);
z0 = zeros(1,Nr+1);
b0 = zeros(1,Nr);

bq1 = zeros(Nt,Nr);
bq2 = zeros(Nt,Nr);
bq3 = zeros(Nt,Nr);
bq4 = zeros(Nt,Nr);
bq5 = zeros(Nt,Nr);
dbdt = zeros(Nt,Nr);
dbdz = zeros(Nt,Nr);

zq1 = zeros(Nt,Nr+1);
zq2 = zeros(Nt,Nr+1);
zq3 = zeros(Nt,Nr+1);
zq4 = zeros(Nt,Nr+1);
dzetadt = zeros(Nt,Nr+1);

d2bdz2 = zeros(1,Nr);

dpsidz = zeros(1,Nr+1);
d2psidz2 = zeros(1,Nr+1);
d2zetadz2 = zeros(1,Nr+1);
dUdz = zeros(1,Nr);
U = zeros(1,Nr);
U_wgrid = zeros(1,Nr+1);

%%% check CLF condition:
CFLz = U0*dt_real/dz_real
CFLx = U0*dt_real/lambda 

%%% Initial condition
% buoy(1,:) = (rand(1,Nr)-1/2)/1e20; %%% Random
% buoy(1,:) = rand(1,Nr)/1e20;       %%% Random 
buoy(1,:) = [1:Nr]/Nr/1e20; %%% Linear 
% buoy(1,:) = 1/1e20;       %%% Constant 
psi(1,:) = 0;
zeta(1,:) = 0;
% zeta(1,:) = (rand(1,Nr)-1/2)/1e30;
% psi(1,:) = (rand(1,Nr)-1/2)/1.79e30;
% for m = 2:Nr-1
%     dpsidz(m)   = (psi(1,m+1)-psi(1,m-1))/2/dz;
%     d2psidz2(m) = (psi(1,m-1)-2*psi(1,m)+psi(1,m+1))/dz^2;
% end
% zeta(1,:) = d2psidz2-kx^2*psi(1,:);

%%% Background tidal velocity (dimensionless)
Atide = zeros(1,Nr);
for m=1:Nr
    if(zz(m)*delta<=Hshear)
        Atide(m) = Uconst + Shear*zz(m)*delta;
    else
        Atide(m) = U0;
    end
end

for m=1:Nr+1
    if(zz_wgrid(m)*delta<=Hshear)
        Atide_wgrid(m) = Uconst + Shear*zz_wgrid(m)*delta;
    else
        Atide_wgrid(m) = U0;
    end
end

% idx_smooth = Nshear-round(Nshear/4):Nshear+4*round(Nshear/4);
idx_smooth = 1:Nr;
% Atide(idx_smooth) = smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(Atide(idx_smooth)))))))))))))))))))))';
Atide(idx_smooth) = smooth(smooth(smooth(Atide(idx_smooth))))';
Atide_wgrid(idx_smooth) = smooth(smooth(smooth(Atide_wgrid(idx_smooth))))';

% figure()
% plot(zz,Atide);grid;hold on;plot(zz_wgrid,Atide_wgrid);

%%
Utide =cos(tt)'.*Atide/U0;

dUtidedz = zeros(Nt,Nr);
for m = 2:Nr-1
    dUtidedz(:,m) = (Utide(:,m+1)-Utide(:,m-1))/2/dz;
end

% h=figure(1);
% set(h,'Visible', FigureIsVisible);clf;
% plot(Atide,zz*delta);
% grid on;grid minor;
% saveas(h,[expdir 'fig1.png'])
% 
% h=figure(2);
% set(h,'Visible', FigureIsVisible);clf;
% plot(diff(Atide)./diff(zz)/delta,0.5*(zz(1:end-1)+zz(2:end))*delta)
% grid on;grid minor;
% saveas(h,[expdir 'fig2.png'])
% 
% h=figure(3);
% set(h,'Visible', FigureIsVisible);clf;
% pcolor(tt/omega/3600,zz*delta,Utide'*U0);shading flat;colormap redblue; colorbar;
% saveas(h,[expdir 'fig3.png'])
% 
% h=figure(4);
% set(h,'Visible', FigureIsVisible);clf;
% pcolor(tt/omega/3600,zz*delta,dUtidedz'*U0/delta);shading flat;colormap redblue; colorbar;
% saveas(h,[expdir 'fig4.png'])

% close all;

%%

%%%% Integrate non-dimensionalized b and zeta with time
%%%% 1st-order and 2nd-order centered difference
%%%% Fourth-order Runge-Kutta or Third-order Adams-Bashforth
%%%% or Euler forward scheme for time advancement

zspan = [0 1];

% %%% Free-slip boundary condition
% zeta(1,1) = 0; zeta(1,Nr+1) = 0; 

for o=1:Nt-1

    if(rem(o,round(Nt/50))==0)
        o/Nt
    end

    t0 = tt(o);
    b0 = buoy(o,:);
    z0 = zeta(o,:);
    
    tendency;
    psi(o,:) = p0;

    if(useRK4)
        RK4;
    elseif(useAB3)
        %%% Third-order Adams-Bashforth method %%%
        if (o <= 2)
            %%% Use RK4 for the first 2 time steps
            RK4; 
            %%% Use Euler-forward for the first 2 time steps
            % buoy(o+1,:) = buoy(o,:) + dbdt(o,:)*dt;
            % zeta(o+1,:) = zeta(o,:) + dzetadt(o,:)*dt;
        else
            buoy(o+1,:) = buoy(o,:) + dt*( (23/12)*dbdt(o,:)    - (16/12)*dbdt(o-1,:)    + (5/12)*dbdt(o-2,:) );
            zeta(o+1,:) = zeta(o,:) + dt*( (23/12)*dzetadt(o,:) - (16/12)*dzetadt(o-1,:) + (5/12)*dzetadt(o-2,:) );
        end
    else
        %%% Euler forward %%%
        buoy(o+1,:) = dbdt(o,:)*dt;
        zeta(o+1,:) = dzetadt(o,:)*dt;
    end

    % %%% Free-slip boundary condition
    zeta(o+1,1) = 0; zeta(o+1,Nr+1) = 0; 

    %%% No-stress (free-slip) b.c. (du/dz = 0) At the ocean bottom
    % zeta(o,1) = delta/Hshear*cos(t0);

    %%% No-stress (free-slip) b.c. (du/dz = 0) At the upper boundary
    % zeta(o,Nr+1) = 0;
    

end

re_psi = real(psi); re_psi(re_psi==0)=NaN;
re_zeta = real(zeta); re_zeta(re_zeta==0)=NaN;
re_buoy = real(buoy); re_buoy(re_buoy==0)=NaN;
re_dbdz = real(dbdz);

re_bq1 = real(bq1);
re_bq2 = real(bq2);
re_bq3 = real(bq3);
re_bq4 = real(bq4);
re_bq5 = real(bq5);

re_zq1 = real(zq1);
re_zq2 = real(zq2);
re_zq3 = real(zq3);
re_zq4 = real(zq4);



%%% Dimensional veriables
zzd = zz*delta;   
zzd_wgrid = zz_wgrid*delta;    
ttd = tt/omega;
U1 = U0;
if(U0==0)
    U1=0.3;
end
re_psid = re_psi*U1/delta;
re_zetad = re_zeta*U1*delta;
re_buoyd = re_buoy*N^2*sind(topo)/omega;

dbuoydz = zeros(Nt,Nr);
for m = 2:Nr-1
    dbuoydz(:,m) = (re_buoyd(:,m+1)-re_buoyd(:,m-1))/2/dz/delta;
end


plot_tidx = 1:10:Nt;
load_colors;

%%
    re_buoy = real(buoy);
    re_zeta = real(zeta);
    re_psi = real(psi);
    
    %%%%% Floquet stability
    oT = round(Ptide*omega/dt);% The time step after one tidal cycle
    % tidx = 10:20;
    zidx=2:Nshear;
    % zidx = 2:Nr
    tidx = 10:20;
    t2 = 3*oT+tidx;
    t1 = 2*oT+tidx;

    % zidx = 2:Nr-1;
    muk_psi = mean(re_psi(oT+tidx,zidx))./mean(re_psi(tidx,zidx));
    muk_zeta = mean(re_zeta(oT+tidx,zidx))./mean(re_zeta(tidx,zidx));
    muk_buoy = mean(re_buoy(oT+tidx,zidx))./mean(re_buoy(tidx,zidx));
    % muk_psi = (re_psi(t2,zidx))./(re_psi(t1,zidx));
    % muk_zeta = (re_zeta(t2,zidx))./(re_zeta(t1,zidx));
    % muk_buoy = (re_buoy(t2,zidx))./(re_buoy(t1,zidx));
    
    muk_max_buoy = max(abs(muk_buoy))
    muk_mean_buoy = mean(abs(muk_buoy));
    muk_rms_buoy = rms(abs(muk_buoy));

    muk_max_zeta = max(abs(muk_zeta))
    muk_mean_zeta = mean(abs(muk_zeta));
    muk_rms_zeta = rms(abs(muk_zeta));

    muk_max_psi = max(abs(muk_psi))
    muk_mean_psi = mean(abs(muk_psi));
    muk_rms_psi = rms(abs(muk_psi));


%%

DIV = 1;
uuu = U1*real((psi(:,2:Nr+1)-psi(:,1:Nr))/dz);
www = real(1i*kx*psi)*U1;

h=figure(5);
set(h,'color','w','Visible', FigureIsVisible,'Position',[67 346 1015 619]);
clf;
subplot(3,2,1)
pcolor(ttd(plot_tidx)/t1hour,zzd_wgrid,re_psid(plot_tidx,:)');shading flat;colorbar;
colormap(redblue)
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');
title('Streamfunction \psi','Fontsize',fontsize+3);
set(gca,'color',gray);
aaa = max(max(abs(re_psid))/DIV);
if(~isnan(aaa))
clim([-1 1]*aaa)
end

subplot(3,2,2)
pcolor(ttd(plot_tidx)/t1hour,zzd_wgrid,re_zetad(plot_tidx,:)');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');
title('Horizontal vorticity perturbation \zeta','Fontsize',fontsize+3);
set(gca,'color',gray);
aaa = max(max(abs(re_zetad))/DIV);
if(~isnan(aaa))
clim([-1 1]*aaa)
end

subplot(3,2,3)
pcolor(ttd(plot_tidx)/t1hour,zzd,re_buoyd(plot_tidx,:)');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');
title('Buoyancy perturbation b^\prime','Fontsize',fontsize+3);
set(gca,'color',gray);
aaa = max(max(abs(re_buoyd))/DIV);
if(~isnan(aaa))
clim([-1 1]*aaa)
end


subplot(3,2,5)
pcolor(ttd(plot_tidx)/t1hour,zzd,uuu(plot_tidx,:)');shading flat;colorbar;
colormap(redblue)
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Horizontal velocity perturbation u^\prime','Fontsize',fontsize+3);
set(gca,'color',gray);
aaa = max(max((abs(uuu))))/DIV;
if(~isnan(aaa)&& aaa~=0)
clim([-1 1]*aaa)
end

subplot(3,2,6)
pcolor(ttd(plot_tidx)/t1hour,zzd_wgrid,www(plot_tidx,:)');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Vertical velocity w','Fontsize',fontsize+3);
set(gca,'color',gray);
if(kx~=0)
    clim([-1 1]*max(max((abs(www))))/DIV)
end

saveas(h,[expdir 'fig5.png'])


decompose;
close all

clear b0 b_wgrid b_2 b_3 b_4 buoy p0 p0_ugrid psi psi0 sol1 solinit ...
    z0 z_2 z_3 z_4 zeta ...
    d2bdz2 d2psidz2 d2zetadz2 dbdt dbdz dpsidz dUtidedz dzetadt ...
    k_1b k_1z k_2b k_2z k_3b k_3z k_4b k_4z h ...
    bq1 bq2 bq3 bq4 bq5 zq1 zq2 zq3 zq4
save(outputname)





