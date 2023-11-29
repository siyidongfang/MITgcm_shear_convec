%%% 
%%% decompose.m
%%%
%%% Calculate the percentage of integrated forcing, 
%%% following Kaiser & Pratt 2022.
%%% Run this script after numerical.m

% Re
% C
% lambda

pposition =  [153 197 1153 319];
blegend = [0.3120 0.1819 0.1425 0.3858];
zlegend =[0.7692 0.1934 0.1168 0.3031];

blegend1 = [0.1351 0.1662 0.1425 0.3858];
zlegend1 = [0.5741 0.1558 0.1168 0.3031];

legend_b = {'$-Ub^\prime_x$',...
    '$-w^\prime \tilde N^2 \cos\theta$',...
    '$-u^\prime \tilde N^2 \sin\theta$',...
    '$-w^\prime B_z$',...
    'Diffusion'};

legend_zeta = {'$-U \zeta^\prime_x$',...
    '$b^\prime_x\cos\theta$',...
    '$-b^\prime_z\cos\theta$',...
    'Diffusion'};

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

% bq_all = bq1_int+bq2_int+bq3_int+bq4_int+bq5_int;
bq_all = 1;

zq1_int = real((sum(zq1(:,zidx),2)))';
zq2_int = real((sum(zq2(:,zidx),2)))';
zq3_int = real((sum(zq3(:,zidx),2)))';
zq4_int = real((sum(zq4(:,zidx),2)))';
% zq_all = zq1_int+zq2_int+zq3_int+zq4_int;
zq_all = 1;

lw = 2;


h=figure(6);
clf;
set(h,'color','w','Position',pposition,'Visible', FigureIsVisible);
subplot(1,2,1)
% plot(ttd(1:Nt-1)/t1hour,bq_all,'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,bq1_int,'LineWidth',lw);
hold on;
plot(ttd(1:Nt-1)/t1hour,bq2_int,'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,bq3_int,'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,bq4_int,'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,bq5_int,'LineWidth',lw,'Color',gray);
set(gca,'Fontsize',fontsize);
xlabel('Time (hours)')
legend(legend_b,'Interpreter','Latex','Fontsize',fontsize+3,'Position',blegend1);
grid on;grid minor;
title('Buoyancy budget')

subplot(1,2,2)
% semilogy(ttd(1:Nt-1)/t1hour,zq_all,'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,(zq1_int),'LineWidth',lw);
hold on;
plot(ttd(1:Nt-1)/t1hour,(zq2_int),'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,(zq3_int),'LineWidth',lw);
plot(ttd(1:Nt-1)/t1hour,(zq4_int),'LineWidth',lw,'Color',gray);
set(gca,'Fontsize',fontsize);xlabel('Time (hours)')
legend(legend_zeta,'Interpreter','Latex','Fontsize',fontsize+3,'Position',zlegend1);
grid on;grid minor;
title('Horizontal vorticity budget')

saveas(h,[expdir 'fig6.png'])


h=figure(7);
clf;
set(h,'color','w','Position',pposition,'Visible', FigureIsVisible);
subplot(1,2,1)
semilogy(ttd(1:Nt-1)/t1hour,abs(bq1_int),'LineWidth',lw);
hold on;
semilogy(ttd(1:Nt-1)/t1hour,abs(bq2_int),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,abs(bq3_int),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,abs(bq4_int),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,abs(bq5_int),'LineWidth',lw,'Color',gray);
set(gca,'Fontsize',fontsize);xlabel('Time (hours)')
legend(legend_b,'Interpreter','Latex','Fontsize',fontsize+3,'Position',blegend);
grid on;grid minor;
title('Buoyancy budget (absolute value, log axis)')


subplot(1,2,2)
semilogy(ttd(1:Nt-1)/t1hour,abs(zq1_int),'LineWidth',lw);
hold on;
semilogy(ttd(1:Nt-1)/t1hour,abs(zq2_int),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,abs(zq3_int),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,abs(zq4_int),'LineWidth',lw,'Color',gray);
set(gca,'Fontsize',fontsize);xlabel('Time (hours)')
legend(legend_zeta,'Interpreter','Latex','Fontsize',fontsize+3,'Position',zlegend);
grid on;grid minor;
title('Horizontal vorticity budget (absolute value, log axis)')

saveas(h,[expdir 'fig7.png'])


wins = round(Nt/NTtide); % window size for move mean

h=figure(8);
clf;
set(h,'color','w','Position',pposition,'Visible', FigureIsVisible);
subplot(1,2,1)
% semilogy(ttd(1:Nt-1)/t1hour,bq_all,'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(bq1_int),wins),'LineWidth',lw);
hold on;
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(bq2_int),wins),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(bq3_int),wins),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(bq4_int),wins),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(bq5_int),wins),'LineWidth',lw,'Color',gray);
set(gca,'Fontsize',fontsize);xlabel('Time (hours)')
legend(legend_b,'Interpreter','Latex','Fontsize',fontsize+3,'Position',blegend);
grid on;grid minor;
title('Buoyancy (absolute value, move mean)')

subplot(1,2,2)
% semilogy(ttd(1:Nt-1)/t1hour,zq_all,'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(zq1_int),wins),'LineWidth',lw);
hold on;
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(zq2_int),wins),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(zq3_int),wins),'LineWidth',lw);
semilogy(ttd(1:Nt-1)/t1hour,movmean(abs(zq4_int),wins),'LineWidth',lw,'Color',gray);
set(gca,'Fontsize',fontsize);xlabel('Time (hours)')
legend(legend_zeta,'Interpreter','Latex','Fontsize',fontsize+3,'Position',zlegend);
grid on;grid minor;
title('Horizontal vorticity (absolute value, move mean)')

saveas(h,[expdir 'fig8.png'])



% bq1_int = real(cumsum(sum(bq1(:,zidx),2)))';
% bq2_int = real(cumsum(sum(bq2(:,zidx),2)))';
% bq3_int = real(cumsum(sum(bq3(:,zidx),2)))';
% bq4_int = real(cumsum(sum(bq4(:,zidx),2)))';
% bq5_int = real(cumsum(sum(bq5(:,zidx),2)))';
% bq_all = 1;
% 
% zq1_int = real(cumsum(sum(zq1(:,zidx),2)))';
% zq2_int = real(cumsum(sum(zq2(:,zidx),2)))';
% zq3_int = real(cumsum(sum(zq3(:,zidx),2)))';
% zq4_int = real(cumsum(sum(zq4(:,zidx),2)))';
% zq_all = 1;
% 
% 
% figure(9)
% clf;set(gcf,'color','w','Position',pposition);
% subplot(1,2,1)
% % semilogy(ttd(1:Nt-1)/t1hour,bq_all,'LineWidth',lw);
% semilogy(ttd(1:Nt-1)/t1hour,abs(bq1_int),'LineWidth',lw);
% hold on;
% semilogy(ttd(1:Nt-1)/t1hour,abs(bq2_int),'LineWidth',lw);
% semilogy(ttd(1:Nt-1)/t1hour,abs(bq3_int),'LineWidth',lw);
% semilogy(ttd(1:Nt-1)/t1hour,abs(bq4_int),'LineWidth',lw);
% semilogy(ttd(1:Nt-1)/t1hour,abs(bq5_int),'LineWidth',lw,'Color',gray);
% set(gca,'Fontsize',fontsize);xlabel('Time (hours)')
% legend(legend_b,'Interpreter','Latex','Fontsize',fontsize+3,'Position',blegend);
% grid on;grid minor;
% 
% subplot(1,2,2)
% % semilogy(ttd(1:Nt-1)/t1hour,zq_all,'LineWidth',lw);
% semilogy(ttd(1:Nt-1)/t1hour,abs(zq1_int),'LineWidth',lw);
% hold on;
% semilogy(ttd(1:Nt-1)/t1hour,abs(zq2_int),'LineWidth',lw);
% semilogy(ttd(1:Nt-1)/t1hour,abs(zq3_int),'LineWidth',lw);
% semilogy(ttd(1:Nt-1)/t1hour,abs(zq4_int),'LineWidth',lw,'Color',gray);
% set(gca,'Fontsize',fontsize);xlabel('Time (hours)')
% legend(legend_zeta,'Interpreter','Latex','Fontsize',fontsize+3,'Position',zlegend);
% grid on;grid minor;

