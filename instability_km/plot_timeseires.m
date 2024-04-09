

% mz_t = m0-rs*st*kx;
% wn2 = kx^2 + mz_t.^2;

% figure(10)
% clf;set(gcf,'Color','w')
% subplot(3,1,1)
% plot(tt/3600/12,mz_t,'LineWidth',2);grid on;grid minor;xlim([0 10]);set(gca,'FontSize',20)
% % xlabel('Time (tidal cycles)')
% title('$m(t)=m_0-\frac{\Lambda}{\omega}\sin(\omega t) k$','Interpreter','latex')
% 
% subplot(3,1,2)
% plot(tt/3600/12,atand(mz_t/kx),'LineWidth',2);grid on;grid minor;xlim([0 10]);set(gca,'FontSize',20)
% % xlabel('Time (tidal cycles)')
% title('$\arctan(\frac{m(t)}{k})$ (degrees)','Interpreter','latex')
% 
% subplot(3,1,3)
% plot(tt/3600/12,wn2,'LineWidth',2);grid on;grid minor;xlim([0 10]);set(gca,'FontSize',20)
% xlabel('Time (tidal cycles)')
% title('$m^2(t)+k^2$','Interpreter','latex')



%%

pe = re_buoy.^2;
ke = 0.5*(re_uuu.^2+re_www.^2);
kew = 0.5*(re_www.^2);
%%% To match Radko (2019) Eq.(19)
pe = pe/4; %%% To match Radko (2019) Eq.(19)
ke = 0.5*((real(-1i*m0*psi)).^2+re_www.^2);
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
clf;set(gcf,'Color','w','Position',[0 1 1467 864])
subplot(2,2,1)
% plot(xxplot*3600/1000,yyplot-yyplot(1),'LineWidth',2)
plot(xxplot/12,yyplot-yyplot(1),'LineWidth',2)
hold on;grid on;grid minor;
% plot(xxplot(fit_span), y_fit(fit_span));
hold off;
% xlabel('Dimensionless $t = t^\star/\tau\ (\tau=10^3\,s)$','Interpreter','latex')
xlabel('Time (tidal cycles)')
set(gca,'FontSize',20)
% xlim([0 2000])
xlim([0 40])
ylabel('$0.5\ln(e)$','Interpreter','latex')
title('Energy')
ylocation = max(yyplot(1:round(Nt/10))-yyplot(1));
text(1,ylocation,['\Lambda=' num2str(shear) ' s^{-1}'],'Color','red','FontSize',50)
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

subplot(2,2,2)
s1_plot = abs(dbdz_fft(1:Nt/2)).^2;
s1_plot = s1_plot/max(s1_plot);

s2_plot = abs(b_fft(1:Nt/2)).^2;
s2_plot = s2_plot/max(s2_plot);

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
subplot(2,2,3)
plot(tt(tplot)/43200,real(buoy(tplot)),'LineWidth',2)
hold on;
plot(tt(tplot)/43200,real(www(tplot))/300,'LineWidth',2)
grid on;grid minor;
xlabel('Time (tidal cycles)')
set(gca,'FontSize',20)
xlim([94 100])
legend('Buoyancy perturbation','Vertical velocity','Position',[0.1395 0.3828 0.1626 0.0583])

% figure(2)
% clf;
% plot(tt(tplot)/43200,pe(tplot));
% hold on;
% plot(tt(tplot)/43200,kew(tplot)/1e5);
% ylim([0 max(pe)/10])


% figure(3)
% clf;set(gcf,'Color','w')
% lb = plot(tt(tplot)/43200,log10_dbdz(tplot),'LineWidth',2);
% hold on;
% lB = plot(tt(tplot)/43200,log10_dBdz(tplot)+log10_dB0dz(tplot),'LineWidth',2);
% ltotal = plot(tt(tplot)/43200,log10_dBdz(tplot)+log10_dbdz(tplot)+log10_dB0dz(tplot),'-.','LineWidth',4);
% lu = plot(tt(tplot)/43200,ct(tplot)/1e6,':','LineWidth',2);
% set(gca,'Fontsize',20);grid on;grid minor;
% xlabel('Time (tidal cycles)')
% axis tight
% % xlim([2.5 5.5])
% legend([lu lb lB ltotal],'Tidal velocity',...
%     'db^\prime/dz','dB_{background}/dz','db_{total}/dz')
% % ylim([-6 6]*1e-6)


%%
period_omega = 2*pi/omega;

lB_plot = dBdz(tplot)*cosd(topo)+dB0dz(tplot)*cosd(topo);
dbdz = dbdz*(max(lB_plot)/max(dbdz));
lB_plot = lB_plot/max(lB_plot);
lb_plot = dbdz(tplot)*cosd(topo)/max(dbdz(tplot)*cosd(topo));
lu_plot = ct(tplot)/max(ct(tplot));
lw_plot = real(www(tplot))/max(real(www(tplot)));
ltotal_plot = dBdz(tplot)*cosd(topo)+dbdz(tplot)*cosd(topo)+dB0dz(tplot)*cosd(topo);
ltotal_plot = ltotal_plot/max(ltotal_plot);

subplot(2,2,4)
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
    'db^\prime/dz','dB_{background}/dz','db_{total}/dz','w',...
    'Position', [0.5756 0.1144 0.1221 0.1574])
title('Normalized values')
% ylim([-6 6]*1e-6)
% end