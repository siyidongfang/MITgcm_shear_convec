
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
R_all = [0:0.05:6];
nR = length(R_all);   
om1 = zeros(1,nR);
om2 = zeros(1,nR);
om3 = zeros(1,nR);
nom1 = zeros(1,nR);
nom2 = zeros(1,nR);
nom3 = zeros(1,nR);
sigma_new_idx = zeros(nR,Nt);
sigma_new = zeros(nR,Nt);
topo=4;

% for o=1:3:nR
for o=100

    R = R_all(o);
    r = cosd(topo)/sind(topo)-N_hat*sqrt(R); %%% m0/k0
    
    A = N_hat^2./(1+N_hat^2*R*(r/N_hat/sqrt(R)-sin(tt_hat)).^2);
    A_harmonic = (1/(sum(1./sqrt(A))/length(A))).^2;

    beta1 = (cosd(topo)-r*sind(topo))^2;
    beta2 = N_hat*sqrt(R)*sind(topo)*(cosd(topo)-r*sind(topo));
    dbeta_dthat = beta2*cos(tt_hat);
    beta = beta1+beta2*sin(tt_hat);
    
    % ana_imag = 0.5*sqrt(4*A.*beta-(cos(tt_hat)./(1+sin(tt_hat))).^2);
    % ana_real = 0.5*(cos(tt_hat)./(1+sin(tt_hat)));
    % ana_imag2 = 0.5*sqrt(4*A.*beta-(dbeta_dthat/beta).^2);
    % ana_real2 = 0.5*(dbeta_dthat./beta);
    % 
    % figure(10)
    % clf;
    % % plot(tt_hat/2/pi,(dbeta_dthat./beta).^2,'g-','LineWidth',2)
    % % hold on;
    % % plot(tt_hat/2/pi,4*A.*beta,'b-','LineWidth',2)
    % plot(tt_hat/2/pi,(dbeta_dthat./beta).^2-4*A.*beta,'g-','LineWidth',2)
    % ylim([-100 100])
    % grid on;grid minor;


    sigma1_imag = imag(0.5*(dbeta_dthat./beta+sqrt((dbeta_dthat/beta)^2-4*A.*beta)));
    sigma2_imag = imag(0.5*(dbeta_dthat./beta-sqrt((dbeta_dthat/beta)^2-4*A.*beta)));
    sigma1_real = real(0.5*(dbeta_dthat./beta+sqrt((dbeta_dthat/beta)^2-4*A.*beta)));
    sigma2_real = real(0.5*(dbeta_dthat./beta-sqrt((dbeta_dthat/beta)^2-4*A.*beta)));

    save('real_component.mat','tt_hat','sigma1_real')

    % figure(11)
    % hold on;
    % plot(tt_hat/2/pi,sigma1_real,'LineWidth',2)
    % ylim([-1 1]*2e3);grid on;grid minor;

    tt_F = 0:dt:2*pi*Ntide*7;
    xi_F = [sigma1_real(1:end-1) sigma1_real(1:end-1) sigma1_real(1:end-1) sigma1_real(1:end-1) sigma1_real(1:end-1) sigma1_real(1:end-1) sigma1_real];
    F = exp(xi_F .* tt_F);

    figure(12)
    hold on;
    plot(tt_F/2/pi,log10(F))
    ylim([-200 200])
    grid on;grid minor;

    sigma = sigma1_imag;
    sigma(sigma==0)=NaN;
    sigma1_imag_harmonic(o) = 1/(sum(1./sigma,'omitnan')/length(sigma));

    for n=1:Nt
        om = sigma1_imag(n);
        if(om>=1)
           om1(o) = om1(o) + 1/om;
           nom1(o) = nom1(o)+1;
           sigma_new_idx(o,n)=1;
        elseif(om<1 && om>=sigma1_imag_harmonic(o))
           om2(o) = om2(o) + 1/om;
           nom2(o) = nom2(o)+1;
           sigma_new_idx(o,n)=2;
        elseif(om>0 && om<sigma1_imag_harmonic(o))
           om3(o) = om3(o) + 1/om;
           nom3(o) = nom3(o)+1;
           sigma_new_idx(o,n)=3;
        elseif(om==0)
           sigma_new_idx(o,n)=3;
        end
    end   

    om1(o) = om1(o)/nom1(o);
    om2(o) = om2(o)/nom2(o);
    om3(o) = om3(o)/nom3(o);
    om1(o) = 1/om1(o);
    om2(o) = 1/om2(o);
    om3(o) = 1/om3(o);

    o1=om1(o);
    o2=om2(o);
    o3=om3(o);

    ntransition(o) = sum(abs(diff(sigma_new_idx(o,:))));
    transition_idx = find(diff(sigma_new_idx(o,:))~=0);

    for n=1:Nt
        aa=sigma_new_idx(o,n);
        if(aa==1)
            sigma_new(o,n)=o1;
        elseif(aa==2)
            sigma_new(o,n)=o2;
        elseif(aa==3)
            sigma_new(o,n)=o3;
        end
    end

    % figure(1)
    % clf;set(gcf,'Color','w','Position',[41 146 831 753])
    % plot(tt_hat/2/pi,sqrt(A),'b-','LineWidth',2)
    % hold on;
    % plot(tt_hat/2/pi,sqrt(A_harmonic)*ones(1,length(tt_hat)),'b:','LineWidth',2)
    % plot(tt_hat/2/pi,sigma1_imag,'r-','LineWidth',2)
    % plot(tt_hat/2/pi,sigma1_imag_harmonic(o)*ones(1,length(tt_hat)),'r:','LineWidth',2)
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
    % hold on
    % plot(tt_hat,sigma2_real,'--')



    % figure(4);clf;set(gcf,'Color','w','Position',[41 146 1275 428])
    % plot(tt_hat/2/pi,sigma1_imag,'LineWidth',2);
    % hold on;
    % plot(tt_hat/2/pi,sigma_new(o,:),'LineWidth',2);
    % plot(tt_hat/2/pi,sigma1_imag_harmonic(o)*ones(1,length(tt_hat)),'--','LineWidth',1.5);
    % plot(tt_hat/2/pi,ones(1,length(tt_hat)),'--','LineWidth',2,'Color','k');
    % grid on;grid minor;set(gca,'Fontsize',20);
    % xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');
    % title(['$\sigma(\hat t)$, R=' num2str(R,2)],'Interpreter','latex');
    % set(gca,'YScale', 'log');
    % legend('$\sigma(\hat t)$','Idealized $\sigma$','harmonic mean of $\sigma(\hat t)$',...
    %     'Interpreter','latex','Position',[0.1447 0.1787 0.1565 0.1647]);
    % ylim([1e-5 10])

end


