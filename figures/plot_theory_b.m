clear all;close all

load('../instability_km/budget/topo4_output_Ri3.mat')

fontsize = 18;

Nall = length(re_uuu);
Nidx = Nall/140*1+1:Nall/140*2;
bb = re_buoy(Nidx);
tt = tt(Nidx);
tt = tt-tt(1);
angle = tt/Ptide*2*pi;

% bb = bb/abs(bb(1));
bb = bb*1e69;
x_plot = bb.*cos(angle);
y_plot = bb.*sin(angle);

yyp =-100:100;
xxp = zeros(1,length(yyp));

figure(1);clf;set(gcf,'Color','w')
plot(xxp,yyp,'k--','LineWidth',0.5)
hold on;
plot(yyp,xxp,'k--','LineWidth',0.5)
s1=scatter(x_plot(1:1080),y_plot(1:1080),40,angle(1:1080)/2/pi*12,'.');
s2=scatter(x_plot(1080:end),y_plot(1080:end),40,angle(1080:end)/2/pi*12,'o');
grid on;grid minor;
mycolormap1=jet;
mycolormap2=sky;
mycolormap2=flip(mycolormap2);
mycolormap=[mycolormap1;mycolormap2];
colormap(mycolormap)
h1=colorbar;
set(get(h1,'Title'),'String','(hour)','interpreter','latex','FontSize',fontsize);
% mycolormap =jet;
% mycolormap=mycolormap(5:end-5,:);
% mycolormap = [mycolormap(1:round(length(mycolormap)/2)-10,:);mycolormap(length(mycolormap):-1:round(length(mycolormap)/2+10),:)];
xlim([-8 4])
ylim([-5 3])
clim([0 12])
pbaspect([12 8 1])
legend([s1 s2],'First half of a period','Second half of a period','Interpreter','latex')
legend('boxoff')
xlabel('$b^\prime/b^\prime_0$','Interpreter','latex')
ylabel('$b^\prime/b^\prime_0$','Interpreter','latex')
set(gca,'Fontsize',fontsize)
title('Normalized buoyancy perturbation','Interpreter','latex','Fontsize',fontsize+5)

print('-dpng','-r300',['fig_supp_new/theory1.png']);

% figure(2)
% plot(angle/2/pi*12,bb)
%%

clear all;

fontsize = 18;


load('../instability_km/budget/flat_output_Ri3.mat')
Nall = length(re_uuu);
Nidx = Nall/140*1+1:Nall/140*2;
bb = re_buoy(Nidx);
tt = tt(Nidx);
tt = tt-tt(1);
angle = tt/Ptide*2*pi;

% bb = bb/abs(bb(1));
bb = bb*1e69;
x_plot = bb.*cos(angle);
y_plot = bb.*sin(angle);


yyp =-100:100;
xxp = zeros(1,length(yyp));

figure(2);clf;set(gcf,'Color','w','Position', [505 102 1224 859])
plot(xxp,yyp,'k--','LineWidth',0.5)
hold on;
plot(yyp,xxp,'k--','LineWidth',0.5)
s1=scatter(x_plot(1:1080),y_plot(1:1080),40,angle(1:1080)/2/pi*12,'.');
s2=scatter(x_plot(1080:end),y_plot(1080:end),40,angle(1080:end)/2/pi*12,'o');
grid on;grid minor;
mycolormap1=jet;
mycolormap2=jet;
mycolormap=[mycolormap1;mycolormap2];
colormap(mycolormap)
h1=colorbar;
set(get(h1,'Title'),'String','Tidal phase (hour)','interpreter','latex','FontSize',fontsize);
legend([s1 s2],'First half of a period','Second half of a period','Interpreter','latex')
legend('boxoff')
% mycolormap =jet;
% mycolormap=mycolormap(5:end-5,:);
% mycolormap = [mycolormap(1:round(length(mycolormap)/2)-10,:);mycolormap(length(mycolormap):-1:round(length(mycolormap)/2+10),:)];
% colormap(mycolormap)
ylim([-2 4])
xlim([-1 30])
clim([0 12])
pbaspect([31 6 1])
xlabel('$b^\prime/b^\prime_0$','Interpreter','latex')
ylabel('$b^\prime/b^\prime_0$','Interpreter','latex')
set(gca,'Fontsize',fontsize)
title('Normalized buoyancy perturbation','Interpreter','latex','Fontsize',fontsize+5)


print('-dpng','-r300',['fig_supp_new/theory2.png']);



%%

clear all;

fontsize = 18;


load('../instability_km/budget/flat_output_nogrow.mat')
Nall = length(re_uuu);
Nidx = Nall/140*1+1:Nall/140*2;
bb = re_buoy(Nidx);
tt = tt(Nidx);
tt = tt-tt(1);
angle = tt/Ptide*2*pi;

bb = bb/abs(bb(1));
% bb = bb*1e70;
x_plot = bb.*cos(angle);
y_plot = bb.*sin(angle);


yyp =-100:100;
xxp = zeros(1,length(yyp));

figure(3);clf;set(gcf,'Color','w')
plot(xxp,yyp,'k--','LineWidth',0.5)
hold on;
plot(yyp,xxp,'k--','LineWidth',0.5)
scatter(x_plot,y_plot,10,angle/2/pi*12,'filled')
grid on;grid minor;
h1=colorbar;
set(get(h1,'Title'),'String','(hour)','interpreter','latex','FontSize',fontsize);
mycolormap1=pink;
mycolormap2=sky;
mycolormap1=mycolormap1(1:end-40,:);
mycolormap2=mycolormap2(41:end,:);
mycolormap2=flip(mycolormap2);
mycolormap=[mycolormap1;mycolormap2];
colormap(mycolormap)

ylim([-1.1 1.1])
xlim([-1.65 1.65])
clim([0 12])
pbaspect([12 8 1])
yticks([-1 0 1])
xticks([-1 0 1])

xlabel('$b^\prime/b^\prime_0$','Interpreter','latex')
ylabel('$b^\prime/b^\prime_0$','Interpreter','latex')
set(gca,'Fontsize',fontsize)
title('Normalized buoyancy perturbation','Interpreter','latex','Fontsize',fontsize+5)


print('-dpng','-r300',['fig_supp_new/theory3.png']);



