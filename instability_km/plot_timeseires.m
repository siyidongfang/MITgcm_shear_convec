

[y_fit,delta_fit] = polyval(pKE,xxplot,S);
figure(20)
clf;
plot(xxplot,yyplot)
hold on;grid on;grid minor;
plot(xxplot(fit_span), y_fit(fit_span));
hold off;

figure(21)
clf;
plot(xxplot,detrend(yyplot))
hold on;grid on;grid minor;
% plot(xxplot(fit_span), y_fit(fit_span));
hold off;


dbdz = 1i*mz*buoy-1i*kx*buoy*rs.*st;
dbdz = real(dbdz);
dbdz = dbdz*1e104;

period = diff(tt([1 end])) ./ (0:1:Nt/2-1);
freq = 2*pi./period;

dbdz_fft=fft(detrend(dbdz))/length(dbdz);
dbdz_fft = dbdz_fft/median(dbdz_fft);

pe_fft = fft(detrend(pe))/length(pe);
pe_fft = pe_fft/median(pe_fft);

kew_fft = fft(detrend(kew))/length(kew);
kew_fft = kew_fft/median(kew_fft);


% st_fft = fft(detrend(s2t))/length(st);
% st_fft = st_fft/median(st_fft);

figure(5)
loglog(freq/omega,abs(dbdz_fft(1:Nt/2)).^2,'LineWidth',2);
% loglog(freq/omega,abs(kew_fft(1:Nt/2)).^2,'LineWidth',2);
% loglog(period/omega,abs(pe_fft(1:Nt/2)).^2,'LineWidth',2);
% semilogx(period/43200,abs(a(1:Nt/2)).^2,'LineWidth',2);
hold on;
grid on;grid minor;

%%
dBdz = -rs*N^2*ss*st;
dB0dz = N^2*cs*ones(1,Nt);

% tplot = Nt/NTtide*5:Nt/NTtide*10;
tplot = 1:Nt;

if(showfig_dbdz)
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

figure(3)
clf;set(gcf,'Color','w')
lb = plot(tt(tplot)/43200,dbdz(tplot),'LineWidth',2);
hold on;
lB = plot(tt(tplot)/43200,dBdz(tplot)+dB0dz(tplot),'LineWidth',2);
ltotal = plot(tt(tplot)/43200,dBdz(tplot)+dbdz(tplot)+dB0dz(tplot),'-.','LineWidth',4);
lu = plot(tt(tplot)/43200,ct(tplot)/1e6,':','LineWidth',2);
set(gca,'Fontsize',20);grid on;grid minor;
xlabel('Time (tidal cycles)')
axis tight
% xlim([2.5 5.5])
legend([lu lb lB ltotal],'Tidal velocity',...
    'db^\prime/dz','dB_{background}/dz','db_{total}/dz')
% ylim([-6 6]*1e-6)
end