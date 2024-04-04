
pe = re_buoy.^2;
ke = 0.5*(re_uuu.^2+re_www.^2);
kew = 0.5*(re_www.^2);
%%% To match Radko (2019) Eq.(19)
pe = pe/4; %%% To match Radko (2019) Eq.(19)
ke = 0.5*((real(-1i*mz*psi)).^2+re_www.^2);
ke = ke/2; 
kew = kew/2;
fit_span = Nt/NTtide*3:Nt;
xxplot = tt/3600;
yyplot = log(pe/median(pe)+ke/median(ke))/2;
% yyplot = log(ke+pe)/2;
[pKE,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
   pKE(1)/3600*1000   
[y_fit,delta_fit] = polyval(pKE,xxplot,S);
figure(20)
clf;set(gcf,'Color','w')
% plot(xxplot*3600/1000,yyplot+25,'LineWidth',2)
plot(xxplot/12,yyplot+25,'LineWidth',2)
hold on;grid on;grid minor;
% plot(xxplot(fit_span), y_fit(fit_span));
hold off;
xlabel('Dimensionless $t = t^\star/\tau\ (\tau=10^3\,s)$','Interpreter','latex')
set(gca,'FontSize',20)
% xlim([0 2000])
ylabel('$0.5\ln(e)$','Interpreter','latex')
%%
figure(21)
clf;
plot(xxplot,detrend(yyplot))
hold on;grid on;grid minor;
% plot(xxplot(fit_span), y_fit(fit_span));
hold off;


dbdz = 1i*mz*buoy-1i*kx*buoy*rs.*st;
dbdz = real(dbdz);
dbdz = dbdz;

period = diff(tt([1 end])) ./ (0:1:Nt/2-1);
freq = 2*pi./period;

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

figure(5)
clf
loglog(freq/omega,abs(dbdz_fft(1:Nt/2)).^2,'LineWidth',2);
% loglog(freq/omega,abs(kew_fft(1:Nt/2)).^2,'LineWidth',2);
% loglog(freq/omega,abs(energy_fft(1:Nt/2)).^2,'LineWidth',2);
% semilogx(period/43200,abs(a(1:Nt/2)).^2,'LineWidth',2);
hold on;
grid on;grid minor;

%%
dBdz = -rs*N^2*ss*st;
dB0dz = N^2*cs*ones(1,Nt);

% tplot = Nt/NTtide*5:Nt/NTtide*10;
tplot = 1:Nt;

% if(showfig_dbdz)
figure(1)
clf
plot(tt(tplot)/43200,real(buoy(tplot)))
hold on;
plot(tt(tplot)/43200,real(www(tplot))/400)

figure(2)
clf;
plot(tt(tplot)/43200,pe(tplot));
hold on;
plot(tt(tplot)/43200,kew(tplot)/1e5);
ylim([0 max(pe)/10])


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
dbdz = dbdz*(N^2)/max(dbdz)*1000;
figure(3)
clf;set(gcf,'Color','w')
lb = plot(tt(tplot)/period_omega,dbdz(tplot),'LineWidth',2);
hold on;
lB = plot(tt(tplot)/period_omega,dBdz(tplot)+dB0dz(tplot),'LineWidth',2);
ltotal = plot(tt(tplot)/period_omega,dBdz(tplot)+dbdz(tplot)+dB0dz(tplot),'-.','LineWidth',4);
lu = plot(tt(tplot)/period_omega,ct(tplot)/4e5,':','LineWidth',2);
lw = plot(tt(tplot)/period_omega,real(www(tplot))*(N^2)/max(www)*1000,'-.','LineWidth',2);
set(gca,'Fontsize',20);grid on;grid minor;
xlabel('Time (tidal cycles)')
axis tight
xlim([54 60])
legend([lu lb lB ltotal lw],'Tidal velocity',...
    'db^\prime/dz','dB_{background}/dz','db_{total}/dz','w')
% ylim([-6 6]*1e-6)
% end