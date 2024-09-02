
clear all;
addpath ../analysis/functions/
load_colors;

N = 1e-3;
omega = 2*pi/43200;
N_hat = N/omega;

Ntide = 1;
dt = 2*pi/10000;
tt_hat = 0:dt:2*pi*Ntide;
Nt = length(tt_hat);

% R_all = [2.1:0.1:6];
R_all = [0:0.1:6];
nR = length(R_all);   
om1 = zeros(1,nR);
om2 = zeros(1,nR);
om3 = zeros(1,nR);
nom1 = zeros(1,nR);
nom2 = zeros(1,nR);
nom3 = zeros(1,nR);
alpha2 = zeros(1,nR);
alpha3 = zeros(1,nR);
NT1 = zeros(1,nR);
NT3 = zeros(1,nR);
NT5 = zeros(1,nR);

topo=4;

for o=1:nR
% for o=10

    R = R_all(o);
    r = cosd(topo)/sind(topo)-N_hat*sqrt(R); %%% m0/k0
    
    A = N_hat^2./(1+N_hat^2*R*(r/N_hat/sqrt(R)-sin(tt_hat)).^2);
    A_harmonic = (1/(sum(1./sqrt(A))/length(A))).^2;

    beta1 = (cosd(topo)-r*sind(topo))^2;
    beta2 = N_hat*sqrt(R)*sind(topo)*(cosd(topo)-r*sind(topo));
    dbeta_dthat = beta2*cos(tt_hat);
    beta = beta1+beta2*sin(tt_hat);
    
    sigma1_imag = imag(0.5*(dbeta_dthat./beta+sqrt((dbeta_dthat/beta)^2-4*A.*beta)));
    sigma2_imag = imag(0.5*(dbeta_dthat./beta-sqrt((dbeta_dthat/beta)^2-4*A.*beta)));
    sigma1_real = real(0.5*(dbeta_dthat./beta+sqrt((dbeta_dthat/beta)^2-4*A.*beta)));
    sigma2_real = real(0.5*(dbeta_dthat./beta-sqrt((dbeta_dthat/beta)^2-4*A.*beta)));

    sigma = sigma1_imag;
    sigma(sigma==0)=NaN;
    sigma1_imag_harmonic(o) = 1/(sum(1./sigma,'omitnan')/length(sigma));

    if(R>2.4)

        for n=1:Nt
            om = sigma1_imag(n);
            if(om~=0)
                if(om>=1)
                    om1(o) = om1(o) + 1/om;
                    nom1(o) = nom1(o)+1;
                    if(n<round(Nt/4))
                        NT1(o) = NT1(o)+1;
                    elseif(n>round(Nt/4) && n<round(Nt/4*3))
                        NT3(o) = NT3(o)+1;
                    else
                        NT5(o) = NT5(o)+1;
                    end
                else
                    if(n<round(Nt/2))
                        om2(o) = om2(o) + 1/om;
                        nom2(o) = nom2(o)+1;
                    else
                        om3(o) = om3(o) + 1/om;
                        nom3(o) = nom3(o)+1;
                    end
                end
            end
        end

        om1(o) = om1(o)/nom1(o);
        om2(o) = om2(o)/nom2(o);
        om3(o) = om3(o)/nom3(o);
        om1(o) = 1/om1(o);
        om2(o) = 1/om2(o);
        om3(o) = 1/om3(o);
        %%% Portion of a period (0~2*pi) with frequency lower than 1
        alpha2(o) = nom2(o)/Nt;
        alpha3(o) = nom3(o)/Nt;

        o1=om1(o);
        o2=om2(o);
        o3=om3(o);
        al2 = alpha2(o);
        al3 = alpha3(o);
        beta = 1-(al2+al3);
        t1 = NT1(o);
        t2 = t1 + nom2(o);
        t3 = t2 + NT3(o);
        t4 = t3 + nom3(o);
        t5 = Nt;
    
        sigma_new = [o1*ones(1,t1) o2*ones(1,t2-t1) o1*ones(1,t3-t2) o3*ones(1,t4-t3) o1*ones(1,t5-t4)];

    else

        for n=1:Nt
            om = sigma1_imag(n);
            if(om~=0)
                if(n<round(Nt/2))
                    om1(o) = om1(o) + 1/om;
                    nom1(o) = nom1(o)+1;
                else
                    om3(o) = om3(o) + 1/om;
                    nom3(o) = nom3(o)+1;
                end
            end
            om2(o) = NaN;
            nom2(o) = 0;
        end

        om1(o) = om1(o)/nom1(o);
        om2(o) = om2(o)/nom2(o);
        om3(o) = om3(o)/nom3(o);
        om1(o) = 1/om1(o);
        om2(o) = 1/om2(o);
        om3(o) = 1/om3(o);
        %%% Portion of a period (0~2*pi) with frequency lower than 1
        alpha2(o) = nom2(o)/Nt;
        alpha3(o) = nom3(o)/Nt;

        o1=om1(o);
        o2=om2(o);
        o3=om3(o);
        al2 = alpha2(o);
        al3 = alpha3(o);
        beta = 1-(al2+al3);
        t1 = NT1(o);
        t2 = t1 + nom2(o);
        t3 = t2 + NT3(o);
        t4 = t3 + nom3(o);
        t5 = Nt;
    
        sigma_new = [o1*ones(1,t4) o3*ones(1,t5-t4)];



    end






    % figure(1)
    % clf;set(gcf,'Color','w','Position',[41 146 831 753])
    % plot(tt_hat/2/pi,sqrt(A),'b-','LineWidth',2)
    % hold on;
    % plot(tt_hat/2/pi,sqrt(A_harmonic)*ones(1,length(tt_hat)),'b:','LineWidth',2)
    % plot(tt_hat/2/pi,sigma1_imag,'r-','LineWidth',2)
    % plot(tt_hat/2/pi,sigma1_imag_harmonic*ones(1,length(tt_hat)),'r:','LineWidth',2)
    % plot(tt_hat/2/pi,ones(1,length(tt_hat)),'--','Color','k','LineWidth',2);
    % set(gca,'YScale','log')
    % % plot(tt_hat,sigma2_imag)
    % grid on;grid minor;set(gca,'Fontsize',20);
    % xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');
    % title(['$\sqrt(A(\hat t))$ and $\sigma(\hat t)$, R=' num2str(R,2)],'Interpreter','latex');
    % legend('$\sqrt(A(\hat t))$','Harmonic mean of $\sqrt(A(\hat t))$'...
    %     ,'$\sigma(\hat t)$','Harmonic mean of $\sigma(\hat t)$','Interpreter','latex',...
    %     'Position', [0.1636 0.1584 0.3105 0.1285]);

    
    % figure(2)
    % clf;set(gcf,'Color','w','Position',[41 146 831 753])
    % plot(tt_hat,sigma1_real,'--')
    % plot(tt_hat,sigma2_real,'--')



    figure(4);clf;set(gcf,'Color','w','Position',[41 146 1275 428])
    plot(tt_hat/2/pi,sigma1_imag,'LineWidth',2);
    hold on;
    plot(tt_hat/2/pi,sigma_new,'LineWidth',2);
    plot(tt_hat/2/pi,sigma1_imag_harmonic(o)*ones(1,length(tt_hat)),'--','LineWidth',1.5);
    plot(tt_hat/2/pi,ones(1,length(tt_hat)),'--','LineWidth',2,'Color','k');
    grid on;grid minor;set(gca,'Fontsize',20);
    xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');
    title(['$\sigma(\hat t)$, R=' num2str(R,2)],'Interpreter','latex');
    set(gca,'YScale', 'log');
    legend('$\sigma(\hat t)$','Idealized $\sigma$','harmonic mean of $\sigma(\hat t)$',...
        'Interpreter','latex','Position',[0.1447 0.1787 0.1565 0.1647]);
    ylim([1e-5 10])

    nt1(o)=t1;
    nt2(o)=t2;
    nt3(o)=t3;
    nt4(o)=t4;
    nt5(o)=t5;

end



close all;
save('freq_topo.mat')

figure(3)
clf;clf;set(gcf,'Color','w','Position', [313 413 1608 465])
subplot(1,2,1)
plot(R_all,om1,'LineWidth',2)
hold on;
plot(R_all,om2,'LineWidth',2,'Color',orange)
plot(R_all,om3,'LineWidth',2,'Color',green)
set(gca,'YScale', 'log');
grid on;grid minor;set(gca,'Fontsize',20);
% ylim([0.1 10])
xlabel('$R={R_i}_\mathrm{min}^{-1}$','Interpreter','latex');
legend('$\omega_1$ = harmonic mean of $\sigma(\hat t)$ higher than 1',...
    '$\omega_2$ = harmonic mean of $\sigma(\hat t)$ lower than 1 ($0<t<\pi$)',...
    '$\omega_3$ = harmonic mean of $\sigma(\hat t)$ lower than 1 ($\pi<t<2\pi$)',...
    'Interpreter','latex');

subplot(1,2,2)
plot(R_all,alpha2,'LineWidth',2,'Color',orange)
hold on;
plot(R_all,alpha3,'LineWidth',2,'Color',green)
plot(R_all,alpha2+alpha3,'k--','LineWidth',2)
grid on;grid minor;set(gca,'Fontsize',20);
xlabel('$R={R_i}_\mathrm{min}^{-1}$','Interpreter','latex');
title('Portion of a period with $\sigma(\hat t)$ lower than 1','Interpreter','latex');
legend('Portion of $\omega_2$','Portion of $\omega_3$','Portion of $\omega_2$ + Portion of $\omega_3$','Interpreter','latex');


