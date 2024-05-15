


outputname = [expdir 'output.mat'];

fontsize = 16;

%%%% Define constants
nu = 2e-4; %%% Kaiser and Pratt 2022: nu=kappa=2e-6
kappa = 2e-4;
Pr = nu/kappa;
t1hour = 3600;
omega = 2*pi/Ptide;
% delta = Hdepth;
% delta = sqrt(2*nu/omega);
% delta = U0/omega;

%%% Model dimension
Lz = Hdepth/delta;     % dimensionless domain height
Nr = round(Lz/dz);
zz = dz/2:dz:(Nr*dz-dz/2);  % Height above topography
zz_wgrid = 0:dz:((Nr)*dz);
dz_real  = dz*delta;

Hshear = 300;
Lshear = Hshear/delta; % dimensionless vertical scale for velocity shear
Nshear = round(Lshear/dz);
Hshear = zz(Nshear)*delta;
U0 = Uconst + Hshear * Shear;
Re = U0*delta/nu;

NTtide = 1.1;
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

zq1 = zeros(Nt,Nr+1);
zq2 = zeros(Nt,Nr+1);
zq3 = zeros(Nt,Nr+1);
zq4 = zeros(Nt,Nr+1);
dzetadt = zeros(Nt,Nr+1);

An = zeros(Nr-2,Nr-2); %%% Matrix
Dn = zeros(1,Nr-2);

dbdz = zeros(1,Nr);
d2bdz2 = zeros(1,Nr);

dpsidz = zeros(1,Nr+1);
d2psidz2 = zeros(1,Nr+1);
d2zetadz2 = zeros(1,Nr+1);
dUdz = zeros(1,Nr);
U = zeros(1,Nr);

% kx = 2*pi/(lambda/delta); 
C = -kx^2*dz^2-2;

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
Atide_wgrid = zeros(1,Nr+1);
for m=1:Nshear
    Atide(m) = Uconst + Shear*zz(m)*delta;
end
for m=1:Nshear+1
    Atide_wgrid(m) = Uconst + Shear*zz_wgrid(m)*delta;
end
for m=Nshear+1:Nr
    Atide(m) = U0;
end
for m=Nshear+2:Nr+1
    Atide_wgrid(m) = U0;
end
% idx_smooth = Nshear-round(Nshear/4):Nshear+4*round(Nshear/4);
idx_smooth = 1:Nr;
% Atide(idx_smooth) = smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(Atide(idx_smooth)))))))))))))))))))))';
Atide(idx_smooth) = smooth(smooth(smooth(Atide(idx_smooth))))';

Utide =cos(tt)'.*Atide/U0;

dUtidedz = zeros(Nt,Nr);
for m = 2:Nr-1
    dUtidedz(:,m)   = (Utide(:,m+1)-Utide(:,m-1))/2/dz;
end

h=figure(1);
set(h,'Visible', FigureIsVisible);clf;
plot(Atide,zz*delta);
grid on;grid minor;
saveas(h,[expdir 'fig1.png'])

h=figure(2);
set(h,'Visible', FigureIsVisible);clf;
plot(diff(Atide)./diff(zz)/delta,0.5*(zz(1:end-1)+zz(2:end))*delta)
grid on;grid minor;
saveas(h,[expdir 'fig2.png'])

h=figure(3);
set(h,'Visible', FigureIsVisible);clf;
pcolor(tt/omega/3600,zz*delta,Utide'*U0);shading flat;colormap redblue; colorbar;
saveas(h,[expdir 'fig3.png'])

h=figure(4);
set(h,'Visible', FigureIsVisible);clf;
pcolor(tt/omega/3600,zz*delta,dUtidedz'*U0/delta);shading flat;colormap redblue; colorbar;
saveas(h,[expdir 'fig4.png'])

% close all;

%%

%%%% Integrate non-dimensionalized b and zeta with time
%%%% 1st-order and 2nd-order centered difference
%%%% Fourth-order Runge-Kutta or Euler forward scheme for time advancement

% buoy(1,:) = rand(1,Nr)/1e20;
% buoy(1,:) = -[1:Nr]/1e20;


for o=1:Nt-1
   
    % o
    if(rem(o,round(Nt/50))==0)
        o/Nt
    end

    t0 = tt(o);
    b0 = buoy(o,:);
    z0 = zeta(o,:);
    
    tendency;
    psi(o,:) = p0;

    if(useRK4)
    %%% Fourth-order Runge-Kutta method %%%
        k_1b = dbdt(o,:);
        k_1z = dzetadt(o,:);
        % Euler forward predictor advancing dt/2:
        b_2 = buoy(o,:)+0.5*dt*k_1b;
        z_2 = zeta(o,:)+0.5*dt*k_1z;
        t0 = tt(o)+dt/2;
        b0 = b_2;
        z0 = z_2;
        tendency;
        k_2b = dbdt(o,:);
        k_2z = dzetadt(o,:);
        % Euler backward corrector advancing dt/2:
        b_3 = buoy(o,:)+0.5*dt*k_2b;
        z_3 = zeta(o,:)+0.5*dt*k_2z;
        t0 = tt(o)+dt/2;
        b0 = b_3;
        z0 = z_3;
        tendency;
        k_3b = dbdt(o,:);
        k_3z = dzetadt(o,:);
        % Mid-point predictor advancing dt:
        b_4 = buoy(o,:)+dt*k_3b;
        z_4 = zeta(o,:)+dt*k_3z;
        t0 = tt(o)+dt;
        b0 = b_4;
        z0 = z_4;
        tendency;
        k_4b = dbdt(o,:);
        k_4z = dzetadt(o,:);
        % Simpson rule corrector advancing dt:
        buoy(o+1,:) = buoy(o,:) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
        zeta(o+1,:) = zeta(o,:) + (1/6)*(k_1z+2*k_2z+2*k_3z+k_4z)*dt;
    elseif(useAB3)
        %%% Third-order Adams-Bashforth method %%%
        if (o <= 2)
            buoy(o+1,:) = buoy(o,:) + dbdt(o,:)*dt;
            zeta(o+1,:) = zeta(o,:) + dzetadt(o,:)*dt;
        else
            buoy(o+1,:) = buoy(o,:) + dt*( (23/12)*dbdt(o,:)    - (16/12)*dbdt(o-1,:)    + (5/12)*dbdt(o-2,:) );
            zeta(o+1,:) = zeta(o,:) + dt*( (23/12)*dzetadt(o,:) - (16/12)*dzetadt(o-1,:) + (5/12)*dzetadt(o-2,:) );
        end
    else
        %%% Euler forward %%%
        buoy(o+1,:) = dbdt(o,:)*dt;
        zeta(o+1,:) = dzetadt(o,:)*dt;
    end

end

re_psi = real(psi); re_psi(re_psi==0)=NaN;
re_zeta = real(zeta); re_zeta(re_zeta==0)=NaN;
re_buoy = real(buoy); re_buoy(re_buoy==0)=NaN;



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

%
plot_tidx = 1:1:Nt;
load_colors;


%%
    re_buoy = real(buoy);
    re_zeta = real(zeta);
    re_psi = real(psi);
    
    %%%%% Floquet stability
    oT = round(Ptide*omega/dt);% The time step after one tidal cycle
    % tidx = 10:20;
    % zidx=2:Nshear;
    tidx = 10;
    zidx = 2:Nr-1;
    % muk_psi = mean(re_psi(oT+tidx,zidx))./mean(re_psi(tidx,zidx));
    % muk_zeta = mean(re_zeta(oT+tidx,zidx))./mean(re_zeta(tidx,zidx));
    % muk_buoy = mean(re_buoy(oT+tidx,zidx))./mean(re_buoy(tidx,zidx));
    muk_psi = (re_psi(oT+tidx,zidx))./(re_psi(tidx,zidx));
    muk_zeta = (re_zeta(oT+tidx,zidx))./(re_zeta(tidx,zidx));
    muk_buoy = (re_buoy(oT+tidx,zidx))./(re_buoy(tidx,zidx));
    
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
% 
% if(ne>=23)
%     DIV = 1;
% elseif(ne>=10)
%     DIV = 1e3;
% else
%     DIV = 1e5;
% end
DIV = 1e5;

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
pcolor(ttd(plot_tidx)/t1hour,zzd_wgrid,log10(abs(re_psid(plot_tidx,:)))');shading flat;colorbar;
colormap(redblue)
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');
title('Streamfunction log10(abs(\psi))','Fontsize',fontsize+3);
set(gca,'color',gray);


subplot(3,2,3)
pcolor(ttd(plot_tidx)/t1hour,zzd_wgrid,re_zetad(plot_tidx,:)');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');
title('Horizontal vorticity perturbation \zeta','Fontsize',fontsize+3);
set(gca,'color',gray);
aaa = max(max(abs(re_zetad))/DIV);
if(~isnan(aaa))
clim([-1 1]*aaa)
end

subplot(3,2,4)
pcolor(ttd(plot_tidx)/t1hour,zzd_wgrid,log10(abs(re_zetad(plot_tidx,:)))');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');
title('Horizontal vorticity perturbation log10(abs(\zeta))','Fontsize',fontsize+3);
set(gca,'color',gray);


subplot(3,2,5)
pcolor(ttd(plot_tidx)/t1hour,zzd,dbuoydz(plot_tidx,:)');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Stratification perturbation','Fontsize',fontsize+3);
set(gca,'color',gray);
aaa = max(max(abs(re_buoyd))/DIV);
if(~isnan(aaa))
clim([-1 1]*aaa)
end

subplot(3,2,6)
pcolor(ttd(plot_tidx)/t1hour,zzd,log10(abs(dbuoydz(plot_tidx,:)))');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Stratification perturbation log10(abs(db/dz))','Fontsize',fontsize+3);
set(gca,'color',gray);

saveas(h,[expdir 'fig5.png'])



%%

uuu = U1*real((psi(:,2:Nr+1)-psi(:,1:Nr))/dz);
% for m=2:Nr-1
%      dpsidz = (psi(:,m+1)-psi(:,m-1))/2/dz;
%      uuu(:,m) = real(dpsidz)*U1;
% end
% uuu(Nr) = uuu(Nr-1);
www = real(1i*kx*psi)*U1;


h=figure(20);
set(h,'color','w', 'Visible',FigureIsVisible,'Position',[137 403 1077 449]);
clf;
subplot(2,2,1)
pcolor(ttd(plot_tidx)/t1hour,zzd,uuu(plot_tidx,:)');shading flat;colorbar;
colormap(redblue)
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');
title('Horizontal velocity perturbation u^\prime','Fontsize',fontsize+3);
set(gca,'color',gray);
aaa = max(max((abs(uuu))))/DIV;
if(~isnan(aaa)&& aaa~=0)
clim([-1 1]*aaa)
end

subplot(2,2,2)
pcolor(ttd(plot_tidx)/t1hour,zzd,log10(abs(uuu(plot_tidx,:)))');shading flat;colorbar;
colormap(redblue)
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');
title('Horizontal velocity perturbation log10(abs(u^\prime))','Fontsize',fontsize+3);
set(gca,'color',gray);

subplot(2,2,3)
pcolor(ttd(plot_tidx)/t1hour,zzd_wgrid,www(plot_tidx,:)');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Vertical velocity w','Fontsize',fontsize+3);
set(gca,'color',gray);
if(kx~=0)
    clim([-1 1]*max(max((abs(www))))/DIV)
end

subplot(2,2,4)
pcolor(ttd(plot_tidx)/t1hour,zzd_wgrid,log10(abs(www(plot_tidx,:)))');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Vertical velocity log10(abs(w))','Fontsize',fontsize+3);
set(gca,'color',gray);

saveas(h,[expdir 'fig9.png'])
% close all

% save(outputname)




