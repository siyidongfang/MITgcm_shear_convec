

% lx = 100;
% lz = 100;
% kx = 2*pi/lx;
% mz = 2*pi/lz;

kx = 10.^([-5:0.5:2]);
mz = 10.^([-5:0.5:2]);
lx = 2*pi./kx;
lz = 2*pi./mz;

Nkx = length(kx);
Nmz = length(mz);

% z=0:0.5:1;
z = 1;
Nz = length(z);

N = 1e-3;
topo = 4;
Ptide = 43200;
omega = 2*pi/Ptide;
NTtide = 10;
dt = 0.01;

Hshear = 250;
Shear = 1.4e-3;
U0 = Hshear * Shear;

nu = 2e-6; 
kappa = 2e-6;

C1 = U0/Hshear/omega;
C2 = N*sind(topo)/omega;
C3 = nu/Hshear^2/omega;
C4 = kappa/Hshear^2/omega;

Lt = NTtide*43200*omega; % dimensionless simulation time
Nt = round(Lt/dt);
tt = dt:dt:Nt*dt;
ttd = tt/omega;

bq1 = zeros(Nkx,Nmz,Nz,Nt);
bq2 = zeros(Nkx,Nmz,Nz,Nt);
bq3 = zeros(Nkx,Nmz,Nz,Nt);
bq4 = zeros(Nkx,Nmz,Nz,Nt);
bq5 = zeros(Nkx,Nmz,Nz,Nt);
dbdt = zeros(Nkx,Nmz,Nz,Nt);

pq1 = zeros(Nkx,Nmz,Nz,Nt);
pq2 = zeros(Nkx,Nmz,Nz,Nt);
pq3 = zeros(Nkx,Nmz,Nz,Nt);
pq4 = zeros(Nkx,Nmz,Nz,Nt);
dpdt = zeros(Nkx,Nmz,Nz,Nt);

buoy = zeros(Nkx,Nmz,Nz,Nt);
psi = zeros(Nkx,Nmz,Nz,Nt);

buoy(:,:,:,1) = 1e-20;

for k = 1:Nkx
    k
    for m = 1:Nmz
        m
        for d = 1:Nz
            for o=1:Nt-1
                
                t0 = tt(o);
                b0 = buoy(k,m,d,o);
                p0 = psi(k,m,d,o);
            
                %%% Fourth-order Runge-Kutta method %%%
                k_1b = dbdt(k,m,d,o);
                k_1p = dpdt(k,m,d,o);
                % Euler forward predictor advancing dt/2:
                b_2 = buoy(k,m,d,o)+0.5*dt*k_1b;
                p_2 = psi(k,m,d,o)+0.5*dt*k_1p;
            
                t0 = tt(o)+dt/2;
                b0 = b_2;
                p0 = p_2;
                tendency_km;
                k_2b = dbdt(k,m,d,o);
                k_2p = dpdt(k,m,d,o);
            
                % Euler backward corrector advancing dt/2:
                b_3 = buoy(k,m,d,o)+0.5*dt*k_2b;
                p_3 = psi(k,m,d,o)+0.5*dt*k_2p;
            
                t0 = tt(o)+dt/2;
                b0 = b_3;
                p0 = p_3;
                tendency_km;
                k_3b = dbdt(k,m,d,o);
                k_3p = dpdt(k,m,d,o);
            
                % Mid-point predictor advancing dt:
                b_4 = buoy(k,m,d,o)+dt*k_3b;
                p_4 = psi(k,m,d,o)+dt*k_3p;
            
                t0 = tt(o)+dt;
                b0 = b_4;
                p0 = p_4;
                tendency_km;
                k_4b = dbdt(k,m,d,o);
                k_4p = dpdt(k,m,d,o);
                % Simpson rule corrector advancing dt:
                buoy(k,m,d,o+1) = buoy(k,m,d,o) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
                psi(k,m,d,o+1) = psi(k,m,d,o) + (1/6)*(k_1p+2*k_2p+2*k_3p+k_4p)*dt;
            
            end
        end
    end
end


re_buoyd = real(buoy)*N^2*U0*sind(topo)/omega;

www = zeros(Nkx,Nmz,Nz,Nt);
for k=1:Nkx
    www(k,:,:,:) = real(1i*kx(k)*psi(k,:,:,:))*U0;
end
uuu = zeros(Nkx,Nmz,Nz,Nt);
for m=1:Nmz
    uuu(:,m,:,:) = real(-1i*mz(m)*psi(:,m,:,:))*U0;
end

ke = 0.5*(uuu.^2+www.^2);
pe = re_buoyd.^2;

%%% Calculate the growth rate for each wave numbers
growth_ke = zeros(Nkx,Nmz,Nz);
growth_pe = zeros(Nkx,Nmz,Nz);
for k=1:Nkx
    for m=1:Nmz
        for d=1:Nz
            growth_ke(k,m,d) =
            growth_pe(k,m,d) =
        end
    end
end



%%
www_plot = squeeze(www(12,10,1,:));
re_buoyd_plot = squeeze(re_buoyd(12,10,1,:));

figure(1)
clf;set(gcf,'Color','w')
subplot(2,1,1)
plot(ttd/3600/12,re_buoyd_plot,'LineWidth',2);
grid on;grid minor;
axis tight;
set(gca,'fontsize',20)
title('Buoyancy perturbation')
xlabel('Time (tidal cycles)');ylabel('(m/s^2)')
subplot(2,1,2)
plot(ttd/3600/12,www_plot,'LineWidth',2);
grid on;grid minor;
axis tight;
set(gca,'fontsize',20)
title('Vertical velocity')
xlabel('Time (tidal cycles)');ylabel('(m/s)')

%%




