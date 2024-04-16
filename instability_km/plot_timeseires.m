


% dt = 600;
% NTtide = 100;
% omega = 2*pi/43200
% Nt = NTtide/omega/dt;
dt_ri = dt/1000;
tt_ri = dt_ri:dt_ri:Nt*dt;
N2 = N^2;
Ri_inverse = (shear*cos(omega*tt_ri)).^2./(N2*cosd(topo) - N2*sind(topo)/omega*shear*sin(omega*tt_ri));
Ri_min = 1/max(Ri_inverse);  

pe = re_buoy.^2;
ke = 0.5*(re_uuu.^2+re_www.^2);
kew = 0.5*(re_www.^2);
%%% To match Radko (2019) Eq.(19)
pe = pe/4; %%% To match Radko (2019) Eq.(19)
% ke = 0.5*((real(-1i*m0*psi)).^2+re_www.^2);
ke = ke/2; 
kew = kew/2;
fit_span = Nt/NTtide*3:Nt;
xxplot = tt/3600;
yyplot = log(pe/median(pe)+ke/median(ke))/2;
% yyplot = log(ke+pe)/2;
[pKE,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
   pKE(1)/3600*1000   
[y_fit,delta_fit] = polyval(pKE,xxplot,S);


figure(1)
clf;set(gcf,'Color','w','Position',[0 1 1467*1.5 864])
subplot(2,3,1)
% plot(xxplot*3600/1000,yyplot-yyplot(1),'LineWidth',2)
plot(xxplot/12,yyplot-yyplot(1),'LineWidth',2)
hold on;grid on;grid minor;
% plot(xxplot(fit_span), y_fit(fit_span));
hold off;
% xlabel('Dimensionless $t = t^\star/\tau\ (\tau=10^3\,s)$','Interpreter','latex')
xlabel('Time (tidal cycles)')
set(gca,'FontSize',20)
% xlim([0 2000])
% xlim([0 40])
ylabel('$0.5\ln(e)$','Interpreter','latex')
title('Energy')
ylocation = max(yyplot(round(Nt/50):round(Nt/20))-yyplot(1));
if(Diffusion)
    ltext = '$\nu=10^{-5}\,m^2/s,\,\kappa=10^{-6}\,m^2/s$';
else
    ltext = '$\nu=\kappa=0$';
end
text(0,ylocation,{['$\Lambda=$' num2str(shear) ' s$^{-1},\, R_{i,\mathrm{min}}=$ ' num2str(Ri_min,3)],...
    ltext,...
    ['$m_0/k_0$=' num2str(m0/kx,2) ', $\arctan(m_0/k_0)$=' num2str(atand(m0/kx),2) '$^\circ$']},...
    'Color','red','FontSize',35,'Interpreter','latex')
%%
% figure(21)
% clf;
% plot(xxplot,detrend(yyplot))
% hold on;grid on;grid minor;
% % plot(xxplot(fit_span), y_fit(fit_span));
% hold off;


dbdz = 1i*m0*buoy-1i*kx*buoy*rs.*st;
dbdz = real(dbdz);

period = diff(tt([1 end])) ./ (0:1:Nt/2-1);
freq = 2*pi./period;

b_fft = fft(detrend(re_buoy))/length(re_buoy);
b_fft = b_fft/median(b_fft);

re_zeta = real(zeta);
zeta_fft = fft(detrend(re_zeta))/length(re_zeta);
zeta_fft = zeta_fft/median(zeta_fft);

dbdz_fft=fft(detrend(dbdz))/length(dbdz);
dbdz_fft = dbdz_fft/median(dbdz_fft);

pe_fft = fft(detrend(pe))/length(pe);
pe_fft = pe_fft/median(pe_fft);

kew_fft = fft(detrend(kew))/length(kew);
kew_fft = kew_fft/median(kew_fft);

energy_fft = fft(detrend((pe/median(pe)+ke/median(ke))/2))/length(pe);
energy_fft = energy_fft/median(energy_fft);

% st_fft = fft(detrend(s2t))/length(st);
% st_fft = st_fft/median(st_fft);

subplot(2,3,2)
s1_plot = abs(dbdz_fft(1:Nt/2)).^2;
s1_plot = s1_plot/max(s1_plot);

s2_plot = abs(b_fft(1:Nt/2)).^2;
s2_plot = s2_plot/max(s2_plot);

% s2_plot = abs(zeta_fft(1:Nt/2)).^2;
% s2_plot = s2_plot/max(s2_plot);

s3_plot = abs(energy_fft(1:Nt/2)).^2;
s3_plot = s3_plot/max(s3_plot);
s2 = loglog(freq/omega,s2_plot,'LineWidth',2);
hold on;
s3 = loglog(freq/omega,s3_plot,'LineWidth',2);
s1 = loglog(freq/omega,s1_plot,'LineWidth',2);
% semilogx(period/43200,abs(a(1:Nt/2)).^2,'LineWidth',2);
grid on;grid minor;
set(gca,'FontSize',20)
title('Normalized spectra')
xlabel('Frequency/\omega')
legend([s1 s2 s3],'db^\prime/dz','b^\prime','energy')

%%
dBdz = -rs*N^2*ss*st;
dB0dz = N^2*cs*ones(1,Nt);

% tplot = Nt/NTtide*5:Nt/NTtide*10;
tplot = 1:Nt;

% if(showfig_dbdz)
subplot(2,3,4)
plot(tt(tplot)/43200,real(buoy(tplot)),'LineWidth',2)
hold on;
plot(tt(tplot)/43200,real(www(tplot))/300,'LineWidth',2)
grid on;grid minor;
xlabel('Time (tidal cycles)')
set(gca,'FontSize',20)
xlim([NTtide-6 NTtide])
legend('Buoyancy perturbation','Vertical velocity','Position', [0.1358 0.3860 0.1022 0.0561])

% figure(2)
% clf;
% plot(tt(tplot)/43200,pe(tplot));
% hold on;
% plot(tt(tplot)/43200,kew(tplot)/1e5);
% ylim([0 max(pe)/10])



period_omega = 2*pi/omega;

lB_plot = dBdz(tplot)/cosd(topo)+dB0dz(tplot)/cosd(topo);

% dbdz = dbdz*(max(lB_plot)/max(dbdz));
% lB_plot = lB_plot/max(lB_plot);

lb_plot = dbdz(tplot)*cosd(topo);
% lb_plot = lb_plot/max(dbdz(tplot)*cosd(topo));
lu_plot = ct(tplot)/max(ct(tplot))*N^2;
lw_plot = real(www(tplot));
lw_plot = lw_plot/max(real(www(tplot)))*N^2;
ltotal_plot = dBdz(tplot)/cosd(topo)+dbdz(tplot)/cosd(topo)+dB0dz(tplot)/cosd(topo);
% ltotal_plot = ltotal_plot/max(ltotal_plot);

subplot(2,3,5)
lb = plot(tt(tplot)/period_omega,lb_plot,'LineWidth',2);
hold on;
lB = plot(tt(tplot)/period_omega,lB_plot,'LineWidth',2);
ltotal = plot(tt(tplot)/period_omega,ltotal_plot,'-.','LineWidth',4);
lu = plot(tt(tplot)/period_omega,lu_plot,':','LineWidth',2);
lw = plot(tt(tplot)/period_omega,lw_plot,'-.','LineWidth',2);
% lw = plot(tt(tplot)/period_omega,real(www(tplot))*(N^2)/max(real(www)),'-.','LineWidth',2);
set(gca,'Fontsize',20);grid on;grid minor;
xlabel('Time (tidal cycles)')
axis tight
xlim([NTtide-6 NTtide])
legend([lu lb lB ltotal lw],'Tidal velocity',...
    'db^\prime/dz','dB_{background}/dz','db_{total}/dz','w','Position', [0.4144 0.1169 0.0768 0.1574])
% title('Normalized values')
% ylim([-6 6]*1e-6)
% end


%%


subplot(2,3,3)
mz_t = m0-rs*st*kx;
wn2 = kx^2 + mz_t.^2;
wn2 = wn2/max(wn2);

plot(tt/3600/12,5*mz_t,'LineWidth',2);grid on;grid minor;xlim([0 10]);set(gca,'FontSize',20)
hold on;
plot(tt/3600/12,atand(mz_t/kx),'LineWidth',2);grid on;grid minor;xlim([0 10]);set(gca,'FontSize',20)
plot(tt/3600/12,50*wn2,'LineWidth',2);grid on;grid minor;xlim([0 10]);set(gca,'FontSize',20)
xlabel('Time (tidal cycles)')

leg1 = '$5\times m(t)=5\times (m_0-\frac{\Lambda}{\omega}\sin(\omega t) k )$';
leg2 = '$\arctan(\frac{m(t)}{k})$ (degrees)';
leg3 = '$50\times(m^2(t)+k^2)/\max(m^2(t)+k^2)$';
legend(leg1,leg2,leg3,'Interpreter','latex','Fontsize',25,'Position',[0.7219 0.3946 0.1863 0.1158])
ylim([-90 90])

