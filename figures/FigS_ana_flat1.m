
clear;
addpath ../transition/
load('../transition/freq_flat_new.mat')


%%
load_colors;

fontsize = 18;

fg1 = figure(1);
clf;
set(gcf,'Color','w','Position', [100 153 1339 800])

tiledlay = tiledlayout(3,3);


o=18;
R = R_all(o);
o1=om1(o);
o2=om2(o);
al = alpha(o);
beta = 1-al;
t1 = round(Nt*beta/4);
t2 = t1 + round(Nt*al/2);
t3 = t2 + round(Nt*beta/2);
t4 = t3 + round(Nt*al/2);
t5 = Nt;

sigma_new = [o1*ones(1,t1) o2*ones(1,t2-t1) o1*ones(1,t3-t2) o2*ones(1,t4-t3) o1*ones(1,t5-t4)];

A = N_hat^2./(1+N_hat^2*R*sin(tt_hat).^2);

xt = 0.78;
yt = 8;
str = '$R=1.7$';

% l21data=sigma_new;
% l21data(sigma_new<1)=NaN;
% l22data=sigma_new;
% l22data(sigma_new>1)=NaN;

nexttile
plot(tt_hat/2/pi,ones(1,length(tt_hat)),'-','LineWidth',1,'Color',black);
hold on;
l1 = plot(tt_hat/2/pi,sqrt(A),'LineWidth',2,'Color',brown);
l2 =plot(tt_hat/2/pi,sigma_new,'LineWidth',2,'Color',pink);
% l21 =plot(tt_hat/2/pi,l21data,'LineWidth',2,'Color',green);
% l22 =plot(tt_hat/2/pi,l22data,'LineWidth',2,'Color',blue);
l3 = plot(tt_hat/2/pi,sigma_harmonic(o)*ones(1,length(tt_hat)),'--','LineWidth',2,'Color',yellow);
set(gca,'YScale', 'log');
legend([l1 l2 l3],'$\sigma(\hat t)$','Idealized $\sigma$ ($\omega_1$ and $\omega_2$)','Harmonic mean of $\sigma(\hat t)$','Interpreter','latex','FontSize',fontsize+4,'Position',[0.0740 0.8323 0.1988 0.1074]);
% legend([l1 l21 l22 l3],'$\sigma(\hat t)$','$\omega_1=$ harmonic mean of $\sigma(\hat t)\geq 1$','$\omega_2=$ harmonic mean of $\sigma(\hat t)<1$','Harmonic mean of $\sigma(\hat t)$','Interpreter','latex','FontSize',fontsize+2);
legend('boxoff')
ylim([0.5 10])
grid on;grid minor;set(gca,'FontSize',fontsize)
title('Time-dependent frequency $\sigma(\hat t)$','Interpreter','latex','FontSize',fontsize+5);
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);




%%
nexttile
hold on;


al = alpha(o);
beta = 1-al;
t1 = round(Nt*beta/4);
t2 = t1 + round(Nt*al/2);
t3 = t2 + round(Nt*beta/2);
t4 = t3 + round(Nt*al/2);
t5 = Nt;

tt1 = tt_hat(1:t1);
tt2 = tt_hat(t1+1:t2);
tt3 = tt_hat(t2+1:t3);
tt4 = tt_hat(t3+1:t4);
tt5 = tt_hat(t4+1:t5);

t1 = t1/Nt*2*pi;
t2 = t2/Nt*2*pi;
t3 = t3/Nt*2*pi;
t4 = t4/Nt*2*pi;
t5 = t5/Nt*2*pi;

nd=1;
zeta = [];
tt = [];

ncycle = 5;
for k=1:ncycle

if(k==1)
   a1=1;
   p1=0;
else
   a1=a5;
   p1=p5;
end

clear p2 a2 p3 a3 p4 a4 p5 a5
syms  p2 a2 p3 a3 p4 a4 p5 a5
eq1 = a1*cos(o1*((k-1)*2*pi+t1)+p1)-a2*cos(o2*((k-1)*2*pi+t1)+p2)==0;
eq2 = o1*a1*sin(o1*((k-1)*2*pi+t1)+p1)-o2*a2*sin(o2*((k-1)*2*pi+t1)+p2)==0;
S = solve(eq1,eq2);
p2 = double(S.p2(nd));
a2 =double(S.a2(nd));

eq3 = a2*cos(o2*((k-1)*2*pi+t2)+p2)-a3*cos(o1*((k-1)*2*pi+t2)+p3)==0;
eq4 = o2*a2*sin(o2*((k-1)*2*pi+t2)+p2)-o1*a3*sin(o1*((k-1)*2*pi+t2)+p3)==0;
S = solve(eq3,eq4);
p3 = double(S.p3(nd));
a3 =double(S.a3(nd));

eq5 = a3*cos(o1*((k-1)*2*pi+t3)+p3)-a4*cos(o2*((k-1)*2*pi+t3)+p4)==0;
eq6 = o1*a3*sin(o1*((k-1)*2*pi+t3)+p3)-o2*a4*sin(o2*((k-1)*2*pi+t3)+p4)==0;
S = solve(eq5,eq6);
p4 = double(S.p4(nd));
a4 =double(S.a4(nd));

eq7 = a4*cos(o2*((k-1)*2*pi+t4)+p4)-a5*cos(o1*((k-1)*2*pi+t4)+p5)==0;
eq8 = o2*a4*sin(o2*((k-1)*2*pi+t4)+p4)-o1*a5*sin(o1*((k-1)*2*pi+t4)+p5)==0;
S = solve(eq7,eq8);
p5 = double(S.p5(nd));
a5 =double(S.a5(nd));

z1 = a1*cos(o1*((k-1)*2*pi+tt1)+p1);
z2 = a2*cos(o2*((k-1)*2*pi+tt2)+p2);
z3 = a3*cos(o1*((k-1)*2*pi+tt3)+p3);
z4 = a4*cos(o2*((k-1)*2*pi+tt4)+p4);
z5 = a5*cos(o1*((k-1)*2*pi+tt5)+p5);

zeta = [zeta z1 z2 z3 z4 z5];

tt = [tt (k-1)*2*pi+tt_hat];



plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',green)
plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',blue)
plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',green)
plot(k-1+tt4/2/pi,z4,'LineWidth',2,'Color',blue)
plot(k-1+tt5/2/pi,z5,'LineWidth',2,'Color',green)

end

grid on;grid minor;set(gca,'FontSize',fontsize)
title('Horizontal vorticity perturbation $\zeta(\hat t)$','Interpreter','latex','FontSize',fontsize+5);
xt = 3.9;
yt = 1.85;
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);
box on;



ke = 0.5*zeta.^2;
fit_span = 1:Nt*ncycle;
xxplot = tt/2/pi*12; %%% in hours
yyplot = log(ke)/2;
[pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
grow(o)=pp(1);
[y_fit,delta_fit] = polyval(pp,xxplot,S);


nexttile
plot(tt/2/pi,log(ke)/2,'LineWidth',2,'Color',darkgray)
hold on;
plot(xxplot(fit_span)/12, y_fit(fit_span),'--','LineWidth',1);
hold off;
grid on;grid minor;set(gca,'FontSize',fontsize)
title('Turbulent kinetic energy $\sim0.5*\zeta^2(\hat t)$','Interpreter','latex','FontSize',fontsize+5);
xt = 3.9;
yt = 10;
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);
legend('$\log$(TKE)/2','Linear fit','Interpreter','latex','FontSize',fontsize+2,'Position',[0.7166 0.8778 0.0914 0.0656]);
legend('boxoff')
ylim([-13 12])



%%
o=25;
R = R_all(o);
o1=om1(o);
o2=om2(o);
al = alpha(o);
beta = 1-al;
t1 = round(Nt*beta/4);
t2 = t1 + round(Nt*al/2);
t3 = t2 + round(Nt*beta/2);
t4 = t3 + round(Nt*al/2);
t5 = Nt;

sigma_new = [o1*ones(1,t1) o2*ones(1,t2-t1) o1*ones(1,t3-t2) o2*ones(1,t4-t3) o1*ones(1,t5-t4)];

A = N_hat^2./(1+N_hat^2*R*sin(tt_hat).^2);

xt = 0.78;
yt = 8;
str = '$R=2.4$';

nexttile
plot(tt_hat/2/pi,ones(1,length(tt_hat)),'-','LineWidth',1,'Color',black);
hold on;
l1 = plot(tt_hat/2/pi,sqrt(A),'LineWidth',2,'Color',brown);
l2 =plot(tt_hat/2/pi,sigma_new,'LineWidth',2,'Color',pink);
l3 = plot(tt_hat/2/pi,sigma_harmonic(o)*ones(1,length(tt_hat)),'--','LineWidth',2,'Color',yellow);
set(gca,'YScale', 'log');
ylim([0.5 10])
grid on;grid minor;set(gca,'FontSize',fontsize)
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);

%%
nexttile
hold on;


al = alpha(o);
beta = 1-al;
t1 = round(Nt*beta/4);
t2 = t1 + round(Nt*al/2);
t3 = t2 + round(Nt*beta/2);
t4 = t3 + round(Nt*al/2);
t5 = Nt;

tt1 = tt_hat(1:t1);
tt2 = tt_hat(t1+1:t2);
tt3 = tt_hat(t2+1:t3);
tt4 = tt_hat(t3+1:t4);
tt5 = tt_hat(t4+1:t5);

t1 = t1/Nt*2*pi;
t2 = t2/Nt*2*pi;
t3 = t3/Nt*2*pi;
t4 = t4/Nt*2*pi;
t5 = t5/Nt*2*pi;

nd=1;
zeta = [];
tt = [];

ncycle = 5;
for k=1:ncycle

if(k==1)
   a1=1;
   p1=0;
else
   a1=a5;
   p1=p5;
end

clear p2 a2 p3 a3 p4 a4 p5 a5
syms  p2 a2 p3 a3 p4 a4 p5 a5
eq1 = a1*cos(o1*((k-1)*2*pi+t1)+p1)-a2*cos(o2*((k-1)*2*pi+t1)+p2)==0;
eq2 = o1*a1*sin(o1*((k-1)*2*pi+t1)+p1)-o2*a2*sin(o2*((k-1)*2*pi+t1)+p2)==0;
S = solve(eq1,eq2);
p2 = double(S.p2(nd));
a2 =double(S.a2(nd));

eq3 = a2*cos(o2*((k-1)*2*pi+t2)+p2)-a3*cos(o1*((k-1)*2*pi+t2)+p3)==0;
eq4 = o2*a2*sin(o2*((k-1)*2*pi+t2)+p2)-o1*a3*sin(o1*((k-1)*2*pi+t2)+p3)==0;
S = solve(eq3,eq4);
p3 = double(S.p3(nd));
a3 =double(S.a3(nd));

eq5 = a3*cos(o1*((k-1)*2*pi+t3)+p3)-a4*cos(o2*((k-1)*2*pi+t3)+p4)==0;
eq6 = o1*a3*sin(o1*((k-1)*2*pi+t3)+p3)-o2*a4*sin(o2*((k-1)*2*pi+t3)+p4)==0;
S = solve(eq5,eq6);
p4 = double(S.p4(nd));
a4 =double(S.a4(nd));

eq7 = a4*cos(o2*((k-1)*2*pi+t4)+p4)-a5*cos(o1*((k-1)*2*pi+t4)+p5)==0;
eq8 = o2*a4*sin(o2*((k-1)*2*pi+t4)+p4)-o1*a5*sin(o1*((k-1)*2*pi+t4)+p5)==0;
S = solve(eq7,eq8);
p5 = double(S.p5(nd));
a5 =double(S.a5(nd));

z1 = a1*cos(o1*((k-1)*2*pi+tt1)+p1);
z2 = a2*cos(o2*((k-1)*2*pi+tt2)+p2);
z3 = a3*cos(o1*((k-1)*2*pi+tt3)+p3);
z4 = a4*cos(o2*((k-1)*2*pi+tt4)+p4);
z5 = a5*cos(o1*((k-1)*2*pi+tt5)+p5);

zeta = [zeta z1 z2 z3 z4 z5];

tt = [tt (k-1)*2*pi+tt_hat];



plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',green)
plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',blue)
plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',green)
plot(k-1+tt4/2/pi,z4,'LineWidth',2,'Color',blue)
plot(k-1+tt5/2/pi,z5,'LineWidth',2,'Color',green)

end

grid on;grid minor;set(gca,'FontSize',fontsize)
xt = 3.9;
yt = 42;
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);
box on;
ylim([-40 50])



ke = 0.5*zeta.^2;
fit_span = 1:Nt*ncycle;
xxplot = tt/2/pi*12; %%% in hours
yyplot = log(ke)/2;
[pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
grow(o)=pp(1);
[y_fit,delta_fit] = polyval(pp,xxplot,S);


nexttile
plot(tt/2/pi,log(ke)/2,'LineWidth',2,'Color',darkgray)
hold on;
plot(xxplot(fit_span)/12, y_fit(fit_span),'--','LineWidth',1);
hold off;
grid on;grid minor;set(gca,'FontSize',fontsize)
xt = 3.9;
yt = 10;
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);

ylim([-13 12])


o=34;
R = R_all(o);
o1=om1(o);
o2=om2(o);
al = alpha(o);
beta = 1-al;
t1 = round(Nt*beta/4);
t2 = t1 + round(Nt*al/2);
t3 = t2 + round(Nt*beta/2);
t4 = t3 + round(Nt*al/2);
t5 = Nt;

sigma_new = [o1*ones(1,t1) o2*ones(1,t2-t1) o1*ones(1,t3-t2) o2*ones(1,t4-t3) o1*ones(1,t5-t4)];

A = N_hat^2./(1+N_hat^2*R*sin(tt_hat).^2);

xt = 0.78;
yt = 8;
str = '$R=3.3$';

nexttile
plot(tt_hat/2/pi,ones(1,length(tt_hat)),'-','LineWidth',1,'Color',black);
hold on;
l1 = plot(tt_hat/2/pi,sqrt(A),'LineWidth',2,'Color',brown);
l2 =plot(tt_hat/2/pi,sigma_new,'LineWidth',2,'Color',pink);
l3 = plot(tt_hat/2/pi,sigma_harmonic(o)*ones(1,length(tt_hat)),'--','LineWidth',2,'Color',yellow);
set(gca,'YScale', 'log');
ylim([0.4 10])
grid on;grid minor;set(gca,'FontSize',fontsize)
xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);




%%
nexttile
hold on;


al = alpha(o);
beta = 1-al;
t1 = round(Nt*beta/4);
t2 = t1 + round(Nt*al/2);
t3 = t2 + round(Nt*beta/2);
t4 = t3 + round(Nt*al/2);
t5 = Nt;

tt1 = tt_hat(1:t1);
tt2 = tt_hat(t1+1:t2);
tt3 = tt_hat(t2+1:t3);
tt4 = tt_hat(t3+1:t4);
tt5 = tt_hat(t4+1:t5);

t1 = t1/Nt*2*pi;
t2 = t2/Nt*2*pi;
t3 = t3/Nt*2*pi;
t4 = t4/Nt*2*pi;
t5 = t5/Nt*2*pi;

nd=1;
zeta = [];
tt = [];

ncycle = 5;
for k=1:ncycle

if(k==1)
   a1=1;
   p1=0;
else
   a1=a5;
   p1=p5;
end

clear p2 a2 p3 a3 p4 a4 p5 a5
syms  p2 a2 p3 a3 p4 a4 p5 a5
eq1 = a1*cos(o1*((k-1)*2*pi+t1)+p1)-a2*cos(o2*((k-1)*2*pi+t1)+p2)==0;
eq2 = o1*a1*sin(o1*((k-1)*2*pi+t1)+p1)-o2*a2*sin(o2*((k-1)*2*pi+t1)+p2)==0;
S = solve(eq1,eq2);
p2 = double(S.p2(nd));
a2 =double(S.a2(nd));

eq3 = a2*cos(o2*((k-1)*2*pi+t2)+p2)-a3*cos(o1*((k-1)*2*pi+t2)+p3)==0;
eq4 = o2*a2*sin(o2*((k-1)*2*pi+t2)+p2)-o1*a3*sin(o1*((k-1)*2*pi+t2)+p3)==0;
S = solve(eq3,eq4);
p3 = double(S.p3(nd));
a3 =double(S.a3(nd));

eq5 = a3*cos(o1*((k-1)*2*pi+t3)+p3)-a4*cos(o2*((k-1)*2*pi+t3)+p4)==0;
eq6 = o1*a3*sin(o1*((k-1)*2*pi+t3)+p3)-o2*a4*sin(o2*((k-1)*2*pi+t3)+p4)==0;
S = solve(eq5,eq6);
p4 = double(S.p4(nd));
a4 =double(S.a4(nd));

eq7 = a4*cos(o2*((k-1)*2*pi+t4)+p4)-a5*cos(o1*((k-1)*2*pi+t4)+p5)==0;
eq8 = o2*a4*sin(o2*((k-1)*2*pi+t4)+p4)-o1*a5*sin(o1*((k-1)*2*pi+t4)+p5)==0;
S = solve(eq7,eq8);
p5 = double(S.p5(nd));
a5 =double(S.a5(nd));

z1 = a1*cos(o1*((k-1)*2*pi+tt1)+p1);
z2 = a2*cos(o2*((k-1)*2*pi+tt2)+p2);
z3 = a3*cos(o1*((k-1)*2*pi+tt3)+p3);
z4 = a4*cos(o2*((k-1)*2*pi+tt4)+p4);
z5 = a5*cos(o1*((k-1)*2*pi+tt5)+p5);

zeta = [zeta z1 z2 z3 z4 z5];

tt = [tt (k-1)*2*pi+tt_hat];



plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',green)
plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',blue)
plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',green)
plot(k-1+tt4/2/pi,z4,'LineWidth',2,'Color',blue)
plot(k-1+tt5/2/pi,z5,'LineWidth',2,'Color',green)

end

grid on;grid minor;set(gca,'FontSize',fontsize)
xt = 3.9;
yt = 8500;
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);
box on;
xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');




ke = 0.5*zeta.^2;
fit_span = 1:Nt*ncycle;
xxplot = tt/2/pi*12; %%% in hours
yyplot = log(ke)/2;
[pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
grow(o)=pp(1);
[y_fit,delta_fit] = polyval(pp,xxplot,S);


nexttile
plot(tt/2/pi,log(ke)/2,'LineWidth',2,'Color',darkgray)
hold on;
plot(xxplot(fit_span)/12, y_fit(fit_span),'--','LineWidth',1);
hold off;
grid on;grid minor;set(gca,'FontSize',fontsize)
xt = 3.9;
yt = 10;
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);
xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');
ylim([-13 12])


tiledlay.TileSpacing = 'compact';
tiledlay.Padding = 'compact';


AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal')



print('-dpng','-r300',['fig_supp_new/figS_ana_flat1_matlab.png']);
