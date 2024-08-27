
clear all;

N = 1e-3;
omega = 2*pi/43200;
N_hat = N/omega;

Ntide = 1;
dt = 2*pi/10000;
tt_hat = 0:dt:2*pi*Ntide;
Nt = length(tt_hat);

Amax = N_hat^2;


R_all = [0:0.1:20];
nr = length(R_all);   
om1 = zeros(1,nr);
om2 = zeros(1,nr);
nom1 = zeros(1,nr);
nom2 = zeros(1,nr);

for o=1:nr
    R = R_all(o);
    
    Amin = N_hat^2/(1+N_hat^2*R);
    A = N_hat^2./(1+N_hat^2*R*sin(tt_hat).^2);
    
    % figure(1);clf;set(gcf,'Color','w','Position',[41 146 1275 428])
    % plot(tt_hat/2/pi,A,'LineWidth',2);
    % hold on;
    % plot(tt_hat/2/pi,ones(1,length(tt_hat)),'--','LineWidth',2);
    % grid on;grid minor;set(gca,'Fontsize',20);
    % xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');
    % title(['$A(\hat t)$, R=' num2str(R,2)],'Interpreter','latex');
    % set(gca,'YScale', 'log');
    % ylim([0.08 100])
    
    for n=1:Nt
        om = sqrt(A(n));
        if(om>=1)
            %%% Harmonic mean of the frequency higher than 1
            om1(o) = om1(o) + 1/om;
            nom1(o) = nom1(o)+1;
        else
            %%% Harmonic mean of the frequency lower than 1
            om2(o) = om2(o) + 1/om;
            nom2(o) = nom2(o)+1;
        end
    end
    
    om1(o) = om1(o)/nom1(o);
    om2(o) = om2(o)/nom2(o);
    om1(o) = 1/om1(o);
    om2(o) = 1/om2(o);
    
    %%% Portion of a period (0~2*pi) with frequency lower than 1
    alpha(o) = nom2(o)/Nt;
end


figure(2)
clf;clf;set(gcf,'Color','w','Position', [313 413 1608 465])
subplot(1,2,1)
plot(R_all,om1,'LineWidth',2)
hold on;
plot(R_all,om2,'LineWidth',2)
set(gca,'YScale', 'log');
grid on;grid minor;set(gca,'Fontsize',20);
ylim([0.1 10])
xlabel('$R={R_i}_\mathrm{min}^{-1}$','Interpreter','latex');
legend('$\omega_1$ = harmonic mean of frequency $\sqrt{A(\hat t)}$ higher than 1',...
    '$\omega_2$ = harmonic mean of frequency $\sqrt{A(\hat t)}$ lower than 1','Interpreter','latex');

subplot(1,2,2)
plot(R_all,alpha,'LineWidth',2)
grid on;grid minor;set(gca,'Fontsize',20);
xlabel('$R={R_i}_\mathrm{min}^{-1}$','Interpreter','latex');
title('Portion of a period with frequency $\sqrt{A(\hat t)}$ lower than 1','Interpreter','latex');



for o=1:nr
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

    A_new = [o1*ones(1,t1) o2*ones(1,t2-t1) o1*ones(1,t3-t2) o2*ones(1,t4-t3) o1*ones(1,t5-t4)];

    A = N_hat^2./(1+N_hat^2*R*sin(tt_hat).^2);

    % figure(3);clf;set(gcf,'Color','w','Position',[41 146 1275 428])
    % plot(tt_hat/2/pi,A,'LineWidth',2);
    % hold on;
    % plot(tt_hat/2/pi,A_new,'LineWidth',2);
    % plot(tt_hat/2/pi,ones(1,length(tt_hat)),'--','LineWidth',2,'Color','k');
    % grid on;grid minor;set(gca,'Fontsize',20);
    % xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');
    % title(['$A(\hat t)$, R=' num2str(R,2)],'Interpreter','latex');
    % set(gca,'YScale', 'log');
    % legend('$A(\hat t)$','Idealized $A$','Interpreter','latex');
    % ylim([0.04 100])
    
    nt1(o)=t1;
    nt2(o)=t2;
    nt3(o)=t3;
    nt4(o)=t4;
    nt5(o)=t5;
end


%%% Save the data
clear A R n o om Amin A_new t1 t2 t3 t4 t5 o1 o2 beta al
save('om_alpha.mat')