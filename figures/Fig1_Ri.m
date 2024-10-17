
clear;close all;
addpath ../observations/
addpath ../analysis/
load('MAVS2_Ri.mat')
time = time/24; % convert into days

load_colors;

fontsize = 17;
gray = [0.7 0.7 0.7];
lightgray = [249 249 249]/255;


meanN2 = mean(n2,'omitnan');
meanShear = mean(shear_int,'omitnan');

time = time-6.75;
XLIM =[0 8];

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1400 700]);


axesposition =[0.048 0.07 0.85 0.5];


ax1 = subplot('position',axesposition);
annotation('textbox',[0.028 0.602 0.15 0.01],'String','h','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
% rec = rectangle('Position',[0 0  2 1.5e-5],'EdgeColor',lightgray,'FaceColor',lightgray,'Curvature',0.1);
% rec.FaceAlpha = 0.5; 
yyaxis left;
plot(time, n2,'Color',gray,'LineWidth',1);
hold on;
% plot(time, meanN2*ones(1,length(time)),'Color',[0 0.4470 0.7410]);
axis tight
plot(time, smooth_n2,'-','LineWidth',1.5,'Color',black);
ax1.YAxis(1).Color = black;

ylabel('$\overline{\partial_{\tilde z}\mathcal B}^{\tilde z}$ (s$^{-2}$)','interpreter','latex');
ylim([-5 10]*1e-6)
xlim(XLIM)

yyaxis right;
plot(time, shear_int,'LineWidth',1.5,'Color',green);
% plot(time, shear_int,'LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
% hold on;
% plot(time, meanShear*ones(1,length(time)),'Color',[0.8500 0.3250 0.0980]);
ylabel('Linear-fit shear $\Lambda(t)$ (s$^{-1}$)','interpreter','latex');
ylim([-2 1.7]*1e-3)
xlim(XLIM)
grid on;grid minor;
ax1.YAxis(2).Color = green;


% xlabel('Dates','interpreter','latex');
set(gca,'FontSize',fontsize,'TickDir', 'in','TickLength',[0.005, 0.005]);
title('Vertical buoyancy gradient, linear-fit shear, and inverse Richardson number of the large-scale flow at MAVS2','interpreter','latex','FontSize',fontsize+5);
hold off;

xticks([0:2:8])
xticklabels({'2021-07-12','2021-07-14','2021-07-16','2021-07-18','2021-07-20'})


ax2 = axes('Position',axesposition, 'Color', 'none');
yyaxis right
plot(time, 1./Ri,':','Color',brown,'LineWidth',1);
hold on;
plot(time, 1./smooth_Ri,'--','Color',blue,'LineWidth',1.5);
axis tight
set(gca,'FontSize',fontsize,'TickDir', 'in','TickLength',[0.005, 0.005]);
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
ylabel('Inverse Richardson number','Color',blue,'Interpreter','latex')
ylim([0 2])
xlim(XLIM)
ax2.YAxis(2).Color = blue;


% print('-dpng','-r300','fig1/fig1_Ri.png');



%%
figure(2)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1400 700]);

ax2 = axes('Position',axesposition, 'Color', 'none');
yyaxis right
axis tight
set(gca,'FontSize',fontsize,'TickDir', 'in','TickLength',[0.005, 0.005]);
ax2.XTick = [];
ylabel('Inverse Richardson number','Color',blue,'Interpreter','latex')
% ylabel('Inverse $R_{i,\mathrm{obs}}$','Color',blue,'Interpreter','latex')
ylim([0 2])
xlim(XLIM)
ax2.YAxis(2).Color = blue;

% print('-dpng','-r300','fig1/fig1_Ri_axes.png');

