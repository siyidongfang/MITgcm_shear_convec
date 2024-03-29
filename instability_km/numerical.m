%%%%% All variables are dimensional variables

clear all;

%%% Define variables
N = 1e-3;
shear = 1.3e-3;
Ptide = 43200;
omega = 2*pi/Ptide;

% kx_all = 10.^([-5:0.5:1.5]);
% mz_all = 10.^([-5:0.5:1.5]);
% lx_all = 2*pi./kx_all;
% lz_all = 2*pi./mz_all;
% nk = length(kx_all);
% nm = length(mz_all);
mz = 1;
% rw_all = 10.^([-2:0.05:1]);
rw_all = 0.1; %%% kx/mz
nr = length(rw_all);

rs = shear/omega; %%% shear over omega
b00 = 1e-7;
b0 = b00*(rand()+rand()*1i);  %%% Initial condition b(t=0)
topo =4;

NTtide = 10;
dt = 30;
Lt = NTtide*43200; 
Nt = Lt/dt;
tt = dt:dt:Nt*dt;

psi = zeros(1,Nt);
zeta = zeros(1,Nt);
buoy = zeros(1,Nt);
dbdt = zeros(1,Nt);
dzetadt = zeros(1,Nt);

cs = cosd(topo);
ss = sind(topo);

%%% Initial condition
buoy(1) = b0;
psi(1) = 0;
zeta(1) = 0;


for j=1:nr
    j
    rw = rw_all(j); %%% ratio of the wavenumbers kx/mz
    kx = mz*rw;

    %%% Start the loop
    for o=1:Nt-1
        %%% Fourth-order Runge-Kutta method %%%
        t0 = tt(o);
        b0 = buoy(o);
        z0 = zeta(o);
        tendency;
        k_1b = dbdt(o);
        k_1z = dzetadt(o);
        % Euler forward predictor advancing dt/2:
        b_2 = buoy(o)+0.5*dt*k_1b;
        z_2 = zeta(o)+0.5*dt*k_1z;
        t0 = tt(o)+dt/2;
        b0 = b_2;
        z0 = z_2;
        tendency;
        k_2b = dbdt(o);
        k_2z = dzetadt(o);
        % Euler backward corrector advancing dt/2:
        b_3 = buoy(o)+0.5*dt*k_2b;
        z_3 = zeta(o)+0.5*dt*k_2z;
        t0 = tt(o)+dt/2;
        b0 = b_3;
        z0 = z_3;
        tendency;
        k_3b = dbdt(o);
        k_3z = dzetadt(o);
        % Mid-point predictor advancing dt:
        b_4 = buoy(o)+dt*k_3b;
        z_4 = zeta(o)+dt*k_3z;
        t0 = tt(o)+dt;
        b0 = b_4;
        z0 = z_4;
        tendency;
        k_4b = dbdt(o);
        k_4z = dzetadt(o);
    
        % Simpson rule corrector advancing dt:
        buoy(o+1) = buoy(o) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
        zeta(o+1) = zeta(o) + (1/6)*(k_1z+2*k_2z+2*k_3z+k_4z)*dt;
    
    end
    
    ct = cos(omega*tt);
    st = sin(omega*tt);
    a1t = -(kx^2+mz^2+2*kx*mz*rs*st+kx^2*rs^2*st.^2);
    a2t = 1i*kx*(2*rs*ss*st-cs)*N^2 - 1i*mz*N^2*ss;
    a3t = 1i*kx*(cs-rs*ss*st) - 1i*mz*ss;
    
    psi = zeta./a1t;
    www = 1i*kx*psi;
    uuu = -1i*mz*psi-1i*kx*psi*rs.*st;
    
    re_buoy = real(buoy);
    re_uuu = real(uuu);
    re_www = real(www);
    pe = re_buoy.^2;
    ke = 0.5*(re_uuu.^2+re_www.^2);
    kew = 0.5*(re_www.^2);
    
    fit_span = Nt/NTtide*2:Nt;
    xxplot = tt/3600;
    yyplot = log(pe/max(pe)+ke/max(ke))/2;
    [pKE,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
    [y_fit,delta_fit] = polyval(pKE,xxplot,S);
    % figure(3)
    % clf;
    % plot(xxplot,yyplot)
    % hold on;
    % plot(xxplot, y_fit)
    growth(j) = pKE(1);

end

if(nr>1)
figure(20)
clf;set(gcf,'Color','w');
plot(log10(rw_all),growth,'LineWidth',2)
grid on;grid minor;
title('Growth rate (1/hour)')
ylabel('(1/hour)')
xlabel('Wavenumber ratio log10(k_x/m_z)')
set(gca,'fontsize',20)
end


figure(1)
clf
plot(tt/43200,real(buoy))
hold on;
plot(tt/43200,real(www)/400)

figure(2)
clf;
plot(tt/43200,pe);
hold on;
plot(tt/43200,kew/1e5);
ylim([0 max(pe)/10])


dbdz = 1i*mz*buoy-1i*kx*buoy*rs.*st;
dbdz = real(dbdz);

dBdz = -rs*N^2*ss*st;
dB0dz = N^2*cs*ones(1,Nt);

%%
figure(3)
clf;set(gcf,'Color','w')
lb = plot(tt/43200,dbdz,'LineWidth',2);
hold on;
lB = plot(tt/43200,dBdz+dB0dz,'LineWidth',2);
ltotal = plot(tt/43200,dBdz+dbdz+dB0dz,'-.','LineWidth',4);
lu = plot(tt/43200,ct/1e6,':','LineWidth',2);
set(gca,'Fontsize',20);grid on;grid minor;
xlabel('Time (tidal cycles)')
axis tight
% xlim([2.5 5.5])
legend([lu lb lB ltotal],'Tidal velocity',...
    'db^\prime/dz','dB_{background}/dz','db_{total}/dz')









