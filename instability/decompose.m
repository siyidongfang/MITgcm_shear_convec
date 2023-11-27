%%% 
%%% decompose.m
%%%
%%% Calculate the percentage of integrated forcing, 
%%% following Kaiser & Pratt 2022.
%%% Run this script after numerical.m

Re
% C

legend_b = ['',...
    '',...
    '',...
    '',...
    '',...
    ''];

legend_zeta = ['',...
    '',...
    '',...
    '',...
    '',...
    '',...
    ''];

%%% To estimate the relative contribution of each term to the instability:
%%% For different wavenumber kx
%%% Vertically integrate the right-hand-side terms of dbdt, dzetadt
%%% Then cumulatively integrate those terms with time

zidx = 1:Nshear;
bq1_int = real((sum(bq1(:,zidx),2)))';
bq2_int = real((sum(bq2(:,zidx),2)))';
bq3_int = real((sum(bq3(:,zidx),2)))';
bq4_int = real((sum(bq4(:,zidx),2)))';
bq5_int = real((sum(bq5(:,zidx),2)))';
bq6_int = real((sum(bq6(:,zidx),2)))';
% bq_all = bq1_int+bq2_int+bq3_int+bq4_int+bq5_int+bq6_int;
bq_all = 1;

zq1_int = real((sum(zq1(:,zidx),2)))';
zq2_int = real((sum(zq2(:,zidx),2)))';
zq3_int = real((sum(zq3(:,zidx),2)))';
zq4_int = real((sum(zq4(:,zidx),2)))';
zq5_int = real((sum(zq5(:,zidx),2)))';
% zq_all = zq1_int+zq2_int+zq3_int+zq4_int+zq5_int+zq6_int+zq7_int;
zq_all = 1;

lw = 1.5;


figure(6)
clf;set(gcf,'color','w');
subplot(1,2,1)
% semilogy(ttd(1:Nt-1)/t1hour,bq_all,'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,(bq1_int./bq_all),'LineWidth',lw);
hold on;
plot(ttd(1:Nt-1)/t1hour,(bq2_int./bq_all),'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,(bq3_int./bq_all),'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,(bq4_int./bq_all),'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,(bq5_int./bq_all),'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,(bq6_int./bq_all),'LineWidth',lw);
legend('bq1','bq2','bq3','bq4','bq5','bq6')
set(gca,'Fontsize',fontsize);
xlabel('Time (hours)')
grid on;grid minor;

subplot(1,2,2)
% semilogy(ttd(1:Nt-1)/t1hour,zq_all,'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,(zq1_int./bq_all),'LineWidth',lw);
hold on;
plot(ttd(1:Nt-1)/t1hour,(zq2_int./bq_all),'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,(zq3_int./bq_all),'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,(zq4_int./bq_all),'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,(zq5_int./bq_all),'LineWidth',lw);
legend('zq1','zq2','zq3','zq4','zq5','zq6','zq7')
set(gca,'Fontsize',fontsize);
xlabel('Time (hours)')
grid on;grid minor;


%%


bq1_int = real(cumsum(sum(bq1(:,zidx),2)))';
bq2_int = real(cumsum(sum(bq2(:,zidx),2)))';
bq3_int = real(cumsum(sum(bq3(:,zidx),2)))';
bq4_int = real(cumsum(sum(bq4(:,zidx),2)))';
bq5_int = real(cumsum(sum(bq5(:,zidx),2)))';
bq6_int = real(cumsum(sum(bq6(:,zidx),2)))';
% bq_all = bq1_int+bq2_int+bq3_int+bq4_int+bq5_int+bq6_int;
bq_all = 1;

zq1_int = real(cumsum(sum(zq1(:,zidx),2)))';
zq2_int = real(cumsum(sum(zq2(:,zidx),2)))';
zq3_int = real(cumsum(sum(zq3(:,zidx),2)))';
zq4_int = real(cumsum(sum(zq4(:,zidx),2)))';
zq5_int = real(cumsum(sum(zq5(:,zidx),2)))';
% zq_all = zq1_int+zq2_int+zq3_int+zq4_int+zq5_int+zq6_int+zq7_int;
zq_all = 1;


figure(8)
clf;set(gcf,'color','w');
subplot(1,2,1)
% semilogy(ttd(1:Nt-1)/t1hour,bq_all,'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,abs(bq1_int./bq_all),'LineWidth',lw);
hold on;
semilogy(ttd(1:Nt-1)/t1hour,abs(bq2_int./bq_all),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,abs(bq3_int./bq_all),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,abs(bq4_int./bq_all),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,abs(bq5_int./bq_all),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,abs(bq6_int./bq_all),'LineWidth',lw);
legend('bq1','bq2','bq3','bq4','bq5','bq6')
set(gca,'Fontsize',fontsize);
xlabel('Time (hours)')
grid on;grid minor;

subplot(1,2,2)
% semilogy(ttd(1:Nt-1)/t1hour,zq_all,'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,abs(zq1_int./bq_all),'LineWidth',lw);
hold on;
semilogy(ttd(1:Nt-1)/t1hour,abs(zq2_int./bq_all),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,abs(zq3_int./bq_all),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,abs(zq4_int./bq_all),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,abs(zq5_int./bq_all),'LineWidth',lw);
legend('zq1','zq2','zq3','zq4','zq5')
set(gca,'Fontsize',fontsize);
xlabel('Time (hours)')
grid on;grid minor;



wins = round(Nt/NTtide)-1; % window size for move mean

figure(10)
clf;set(gcf,'color','w');
subplot(1,2,1)
% semilogy(ttd(1:Nt-1)/t1hour,bq_all,'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(bq1_int./bq_all),wins),'LineWidth',lw);
hold on;
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(bq2_int./bq_all),wins),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(bq3_int./bq_all),wins),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(bq4_int./bq_all),wins),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(bq5_int./bq_all),wins),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(bq6_int./bq_all),wins),'LineWidth',lw);
legend('bq1','bq2','bq3','bq4','bq5','bq6')
set(gca,'Fontsize',fontsize);
xlabel('Time (hours)')
grid on;grid minor;

subplot(1,2,2)
% semilogy(ttd(1:Nt-1)/t1hour,zq_all,'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(zq1_int./bq_all),wins),'LineWidth',lw);
hold on;
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(zq2_int./bq_all),wins),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(zq3_int./bq_all),wins),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(zq4_int./bq_all),wins),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(zq5_int./bq_all),wins),'LineWidth',lw);
legend('zq1','zq2','zq3','zq4','zq5')
set(gca,'Fontsize',fontsize);
xlabel('Time (hours)')
grid on;grid minor;
