%%%%% All variables are dimensional variables

clear all;
close all;

showfigure = false;
N = sqrt(10)*1e-3;
shear_all = [0:0.1:3]*1e-3;
ns = length(shear_all);
Ptide = 43200;
% omega = 2*pi/Ptide;
omega = 0.1*1e-3;

% mz_all = [0:0.03:6];
% kx_all = [-0.1:0.001:0.1];

% rw_all = 10.^([-2:0.1:-1.2 -1.15:0.05:0.6 1 2 3 4]); %%% kx/mz
% rw_all = 10.^([-3:0.005:1]); %%% kx/mz
rw_all = 10.^([-3:0.02:1]); %%% kx/mz
% nr = length(rw_all);
topo_all = [0:9];

b00 = 1e-70;
b0 = b00*(rand()+rand()*1i);  %%% Initial condition b(t=0)

NTtide = 100;
dt = 600;
Lt = NTtide*43200; 
Nt = Lt/dt;
tt = dt:dt:Nt*dt;

kappa = 1e-7;
nu = 1e-6;

psi = zeros(1,Nt);
zeta = zeros(1,Nt);
buoy = zeros(1,Nt);
dbdt = zeros(1,Nt);
dzetadt = zeros(1,Nt);

%%% Initial condition
buoy(1) = b0;
psi(1) = 0;
zeta(1) = 0;


% for m =1:length(mz_all)
%     m
    % mz = mz_all(m);
    mz = 0.01;

for topo = topo_all(1)

    cs = cosd(topo);
    ss = sind(topo);

% for i=1:ns
for i = 11
    % shear = shear_all(i);
    shear = sqrt(N^2);
    rs = shear/omega; %%% shear over omega

    % for j=1:nr
    % for j=1:length(rw_all)
        % j
        % rw = rw_all(j); %%% ratio of the wavenumbers kx/mz
        rw = 0.015;
        kx = mz*rw;
        % kx = kx_all(j);
    
        %%% Start the loop
        for o=1:Nt-1
            %%% Fourth-order Runge-Kutta method %%%
            t0 = tt(o);
            b0 = buoy(o);
            z0 = zeta(o);
            % tendency_NOdiff;
            tendency;
            k_1b = dbdt(o);
            k_1z = dzetadt(o);
            % Euler forward predictor advancing dt/2:
            b_2 = buoy(o)+0.5*dt*k_1b;
            z_2 = zeta(o)+0.5*dt*k_1z;
            t0 = tt(o)+dt/2;
            b0 = b_2;
            z0 = z_2;
            % tendency_NOdiff;
            tendency;
            k_2b = dbdt(o);
            k_2z = dzetadt(o);
            % Euler backward corrector advancing dt/2:
            b_3 = buoy(o)+0.5*dt*k_2b;
            z_3 = zeta(o)+0.5*dt*k_2z;
            t0 = tt(o)+dt/2;
            b0 = b_3;
            z0 = z_3;
            % tendency_NOdiff;
            tendency;
            k_3b = dbdt(o);
            k_3z = dzetadt(o);
            % Mid-point predictor advancing dt:
            b_4 = buoy(o)+dt*k_3b;
            z_4 = zeta(o)+dt*k_3z;
            t0 = tt(o)+dt;
            b0 = b_4;
            z0 = z_4;
            % tendency_NOdiff;
            tendency;
            k_4b = dbdt(o);
            k_4z = dzetadt(o);
        
            % Simpson rule corrector advancing dt:
            buoy(o+1) = buoy(o) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
            zeta(o+1) = zeta(o) + (1/6)*(k_1z+2*k_2z+2*k_3z+k_4z)*dt;

        end
        
        ct = cos(omega*tt);
        st = sin(omega*tt);
    
        a1t = -(kx^2+mz^2+kx^2*rs^2*st.^2)+2*kx*mz*rs*st;
        
        psi = zeta./a1t;
        www = 1i*kx*psi;
        uuu = -1i*mz*psi-1i*kx*psi*rs.*st;
        
        re_buoy = real(buoy);
        re_uuu = real(uuu);
        re_www = real(www);
        pe = re_buoy.^2;
        ke = 0.5*(re_uuu.^2+re_www.^2);
        kew = 0.5*(re_www.^2);
        % %%% To match Radko (2019) Eq.(19)
        % pe = pe/4; %%% To match Radko (2019) Eq.(19)
        % ke = 0.5*((real(-1i*mz*psi)).^2+re_www.^2);
        % ke = ke/2; 
        % kew = kew/2;
        
        fit_span = Nt/NTtide*3:Nt;
        xxplot = tt/3600;
        yyplot = log(pe/median(pe)+ke/median(ke))/2;
        % yyplot = log(pe+ke)/2;
        [pKE,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
        % [y_fit,delta_fit] = polyval(pKE,xxplot,S);
        % figure(20)
        % clf;
        % plot(xxplot,yyplot)
        % hold on;grid on;grid minor;
        % plot(xxplot(fit_span), y_fit(fit_span));
        % hold off;
        % growth(j) = pKE(1)
        save(['output_Ri1/topo' num2str(topo) '/mz' num2str(mz) 'kx' num2str(kx) 'shear' num2str(shear) '_NOdiff.mat'])
    % end
end

% save(['output_Ri1/growth_topo' num2str(topo) '_mz' num2str(mz) 'kx' num2str(kx) '_NOdiff.mat'])

end

% if(nr>1)
% figure(21);
% % clf
% set(gcf,'Color','w');
% plot(log10(rw_all),growth,'LineWidth',2)
% hold on;
% grid on;grid minor;
% title('Growth rate (1/hour)')
% ylabel('(1/hour)')
% xlabel('Wavenumber ratio log10(k_x/m_z)')
% set(gca,'fontsize',20)
% end


%%%
% if(nr==1)
% plot_timeseires
% end


% end


% save(['output_Ri1/growth_topo' num2str(topo) '_NOdiff_test.mat'])






