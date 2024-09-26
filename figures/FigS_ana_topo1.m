
clear;close all;
addpath ../transition/
load('../transition/freq_topo_3om.mat')
load_colors;

fontsize = 16;

fg1 = figure(1);
clf;
set(gcf,'Color','w','Position', [100 153 1300 1000])

tiledlay = tiledlayout(4,3);


xt = 0.78;
yt = 4;
str = '$R=0.5$';
o=21;

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


nexttile
plot(tt_hat/2/pi,ones(1,length(tt_hat)),'-','LineWidth',1,'Color',black);
hold on;
l1 = plot(tt_hat/2/pi,sigma1_imag,'LineWidth',2,'Color',brown);
l2 = plot(tt_hat/2/pi,sigma_new(o,:),'LineWidth',2,'Color',pink);
l3 = plot(tt_hat/2/pi,sigma1_imag_harmonic(o)*ones(1,length(tt_hat)),'--','LineWidth',2,'Color',yellow);
set(gca,'YScale', 'log');
legend([l1 l2 l3],'$\sigma(\hat t)$','Idealized $\sigma$ ($\omega_1$, $\omega_2$ and $\omega_3$)','$\overline{\sigma}^\mathrm{hm}$, harmonic mean of $\sigma(\hat t)$','Interpreter','latex','FontSize',fontsize+4,'Position', [0.0644 0.7776 0.1890 0.0805]);
legend('boxoff')
ylim([1e-4 10])
grid on;grid minor;set(gca,'FontSize',fontsize)
title('Time-dependent frequency $\sigma(\hat t)$','Interpreter','latex','FontSize',fontsize+5);
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);




%%
nexttile
hold on;

    R = R_all(o);

    o1=om1(o);
    o2=om2(o);
    o3=om3(o);


    ntransition(o);
    transition_idx = find(diff(sigma_new_idx(o,:))~=0);

    if(ntransition(o)==2)
        t1 = transition_idx(1);
        t2 = transition_idx(2);
        t3 = Nt;
        tt1 = tt_hat(1:t1);
        tt2 = tt_hat(t1+1:t2);
        tt3 = tt_hat(t2+1:t3);
        t1 = t1/Nt*2*pi;
        t2 = t2/Nt*2*pi;
        t3 = t3/Nt*2*pi;
    elseif(ntransition(o)==4)
        t1 = transition_idx(1);
        t2 = transition_idx(2);
        t3 = transition_idx(3);
        t4 = transition_idx(4);
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
    elseif(ntransition(o)==6)
        t1 = transition_idx(1);
        t2 = transition_idx(2);
        t3 = transition_idx(3);
        t4 = transition_idx(4);
        t5 = transition_idx(5);
        t6 = transition_idx(6);
        t7 = Nt;
        tt1 = tt_hat(1:t1);
        tt2 = tt_hat(t1+1:t2);
        tt3 = tt_hat(t2+1:t3);
        tt4 = tt_hat(t3+1:t4);
        tt5 = tt_hat(t4+1:t5);
        tt6 = tt_hat(t5+1:t6);
        tt7 = tt_hat(t6+1:t7);
        t1 = t1/Nt*2*pi;
        t2 = t2/Nt*2*pi;
        t3 = t3/Nt*2*pi;
        t4 = t4/Nt*2*pi;
        t5 = t5/Nt*2*pi;
        t6 = t6/Nt*2*pi;
        t7 = t7/Nt*2*pi;
    end

    nd=1;
    zeta = [];
    tt = [];
    
    ncycle = 10;
    for k=1:ncycle

        if(ntransition(o)==2)
            if(k==1)
               a1=1;
               p1=0;
            else
               a1=a3;
               p1=p3;
            end
            
            clear p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7
            syms  p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7

            eq1 = a1*cos(o2*((k-1)*2*pi+t1)+p1)-a2*cos(o3*((k-1)*2*pi+t1)+p2)==0;
            eq2 = o2*a1*sin(o2*((k-1)*2*pi+t1)+p1)-o3*a2*sin(o3*((k-1)*2*pi+t1)+p2)==0;
            S = solve(eq1,eq2);
            p2 = double(S.p2(nd));
            a2 =double(S.a2(nd));
            
            eq3 = a2*cos(o3*((k-1)*2*pi+t2)+p2)-a3*cos(o2*((k-1)*2*pi+t2)+p3)==0;
            eq4 = o3*a2*sin(o3*((k-1)*2*pi+t2)+p2)-o2*a3*sin(o2*((k-1)*2*pi+t2)+p3)==0;
            S = solve(eq3,eq4);
            p3 = double(S.p3(nd));
            a3 =double(S.a3(nd));

            z1 = a1*cos(o2*((k-1)*2*pi+tt1)+p1);
            z2 = a2*cos(o3*((k-1)*2*pi+tt2)+p2);
            z3 = a3*cos(o2*((k-1)*2*pi+tt3)+p3);
            z4=[];
            z5=[];
            z6=[];
            z7=[];

            zeta = [zeta z1 z2 z3 z4 z5 z6 z7];
            
            tt = [tt (k-1)*2*pi+tt_hat];
            plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',blue)
            plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',orange)
            plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',blue)
            
        elseif(ntransition(o)==4)

            if(k==1)
               a1=1;
               p1=0;
            else
               a1=a5;
               p1=p5;
            end
            
            clear p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7
            syms  p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7

            eq1 = a1*cos(o2*((k-1)*2*pi+t1)+p1)-a2*cos(o1*((k-1)*2*pi+t1)+p2)==0;
            eq2 = o2*a1*sin(o2*((k-1)*2*pi+t1)+p1)-o1*a2*sin(o1*((k-1)*2*pi+t1)+p2)==0;
            S = solve(eq1,eq2);
            p2 = double(S.p2(nd));
            a2 =double(S.a2(nd));
            
            eq3 = a2*cos(o1*((k-1)*2*pi+t2)+p2)-a3*cos(o2*((k-1)*2*pi+t2)+p3)==0;
            eq4 = o1*a2*sin(o1*((k-1)*2*pi+t2)+p2)-o2*a3*sin(o2*((k-1)*2*pi+t2)+p3)==0;
            S = solve(eq3,eq4);
            p3 = double(S.p3(nd));
            a3 =double(S.a3(nd));
            
            eq5 = a3*cos(o2*((k-1)*2*pi+t3)+p3)-a4*cos(o3*((k-1)*2*pi+t3)+p4)==0;
            eq6 = o2*a3*sin(o2*((k-1)*2*pi+t3)+p3)-o3*a4*sin(o3*((k-1)*2*pi+t3)+p4)==0;
            S = solve(eq5,eq6);
            p4 = double(S.p4(nd));
            a4 =double(S.a4(nd));
            
            eq7 = a4*cos(o3*((k-1)*2*pi+t4)+p4)-a5*cos(o2*((k-1)*2*pi+t4)+p5)==0;
            eq8 = o3*a4*sin(o3*((k-1)*2*pi+t4)+p4)-o2*a5*sin(o2*((k-1)*2*pi+t4)+p5)==0;
            S = solve(eq7,eq8);
            p5 = double(S.p5(nd));
            a5 =double(S.a5(nd));

            z1 = a1*cos(o2*((k-1)*2*pi+tt1)+p1);
            z2 = a2*cos(o1*((k-1)*2*pi+tt2)+p2);
            z3 = a3*cos(o2*((k-1)*2*pi+tt3)+p3);
            z4 = a4*cos(o3*((k-1)*2*pi+tt4)+p4);
            z5 = a5*cos(o2*((k-1)*2*pi+tt5)+p5);
            z6=[];
            z7=[];

            zeta = [zeta z1 z2 z3 z4 z5 z6 z7];
            
            tt = [tt (k-1)*2*pi+tt_hat];
            plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',green)
            plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',blue)
            plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',orange)
            plot(k-1+tt4/2/pi,z4,'LineWidth',2,'Color',blue)
            plot(k-1+tt5/2/pi,z5,'LineWidth',2,'Color',green)

        elseif(ntransition(o)==6)
            if(k==1)
               a1=1;
               p1=0;
            else
               a1=a7;
               p1=p7;
            end
            
            clear p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7
            syms  p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7

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
            
            eq7 = a4*cos(o2*((k-1)*2*pi+t4)+p4)-a5*cos(o3*((k-1)*2*pi+t4)+p5)==0;
            eq8 = o2*a4*sin(o2*((k-1)*2*pi+t4)+p4)-o3*a5*sin(o3*((k-1)*2*pi+t4)+p5)==0;
            S = solve(eq7,eq8);
            p5 = double(S.p5(nd));
            a5 =double(S.a5(nd));

            eq9 = a5*cos(o3*((k-1)*2*pi+t5)+p5)-a6*cos(o2*((k-1)*2*pi+t5)+p6)==0;
            eq10 = o3*a5*sin(o3*((k-1)*2*pi+t5)+p5)-o2*a6*sin(o2*((k-1)*2*pi+t5)+p6)==0;
            S = solve(eq9,eq10);
            p6 = double(S.p6(nd));
            a6 =double(S.a6(nd));
            
            eq11 = a6*cos(o2*((k-1)*2*pi+t6)+p6)-a7*cos(o1*((k-1)*2*pi+t6)+p7)==0;
            eq12 = o2*a6*sin(o2*((k-1)*2*pi+t6)+p6)-o1*a7*sin(o1*((k-1)*2*pi+t6)+p7)==0;
            S = solve(eq11,eq12);
            p7 = double(S.p7(nd));
            a7 =double(S.a7(nd));

            z1 = a1*cos(o1*((k-1)*2*pi+tt1)+p1);
            z2 = a2*cos(o2*((k-1)*2*pi+tt2)+p2);
            z3 = a3*cos(o1*((k-1)*2*pi+tt3)+p3);
            z4 = a4*cos(o2*((k-1)*2*pi+tt4)+p4);
            z5 = a5*cos(o3*((k-1)*2*pi+tt5)+p5);
            z6 = a6*cos(o2*((k-1)*2*pi+tt6)+p6);
            z7 = a7*cos(o1*((k-1)*2*pi+tt7)+p7);

            zeta = [zeta z1 z2 z3 z4 z5 z6 z7];
            
            tt = [tt (k-1)*2*pi+tt_hat];
            plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',green)
            plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',blue)
            plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',green)
            plot(k-1+tt4/2/pi,z4,'LineWidth',2,'Color',blue)
            plot(k-1+tt5/2/pi,z5,'LineWidth',2,'Color',orange)
            plot(k-1+tt6/2/pi,z6,'LineWidth',2,'Color',blue)
            plot(k-1+tt7/2/pi,z7,'LineWidth',2,'Color',green)

        end
        

    end


grid on;grid minor;set(gca,'FontSize',fontsize)
title('Horizontal vorticity perturbation $\zeta(\hat t)$','Interpreter','latex','FontSize',fontsize+5);
xt = 7.9;
yt = 0.84;
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
xt = 7.9;
yt = 11;
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);
legend('$\log$(TKE)/2','Linear fit','Interpreter','latex','FontSize',fontsize+2,'Position',[0.7166 0.8778 0.0914 0.0656]);
legend('boxoff')
ylim([-14 13])



%%


xt = 0.78;
yt = 4;
str = '$R=1.1$';

o=45;

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


nexttile
plot(tt_hat/2/pi,ones(1,length(tt_hat)),'-','LineWidth',1,'Color',black);
hold on;
l1 = plot(tt_hat/2/pi,sigma1_imag,'LineWidth',2,'Color',brown);
l2 = plot(tt_hat/2/pi,sigma_new(o,:),'LineWidth',2,'Color',pink);
l3 = plot(tt_hat/2/pi,sigma1_imag_harmonic(o)*ones(1,length(tt_hat)),'--','LineWidth',2,'Color',yellow);
set(gca,'YScale', 'log');
ylim([1e-4 12])
grid on;grid minor;set(gca,'FontSize',fontsize)
title('Time-dependent frequency $\sigma(\hat t)$','Interpreter','latex','FontSize',fontsize+5);
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);


%%
nexttile
hold on;

    R = R_all(o);

    o1=om1(o);
    o2=om2(o);
    o3=om3(o);


    ntransition(o);
    transition_idx = find(diff(sigma_new_idx(o,:))~=0);

    if(ntransition(o)==2)
        t1 = transition_idx(1);
        t2 = transition_idx(2);
        t3 = Nt;
        tt1 = tt_hat(1:t1);
        tt2 = tt_hat(t1+1:t2);
        tt3 = tt_hat(t2+1:t3);
        t1 = t1/Nt*2*pi;
        t2 = t2/Nt*2*pi;
        t3 = t3/Nt*2*pi;
    elseif(ntransition(o)==4)
        t1 = transition_idx(1);
        t2 = transition_idx(2);
        t3 = transition_idx(3);
        t4 = transition_idx(4);
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
    elseif(ntransition(o)==6)
        t1 = transition_idx(1);
        t2 = transition_idx(2);
        t3 = transition_idx(3);
        t4 = transition_idx(4);
        t5 = transition_idx(5);
        t6 = transition_idx(6);
        t7 = Nt;
        tt1 = tt_hat(1:t1);
        tt2 = tt_hat(t1+1:t2);
        tt3 = tt_hat(t2+1:t3);
        tt4 = tt_hat(t3+1:t4);
        tt5 = tt_hat(t4+1:t5);
        tt6 = tt_hat(t5+1:t6);
        tt7 = tt_hat(t6+1:t7);
        t1 = t1/Nt*2*pi;
        t2 = t2/Nt*2*pi;
        t3 = t3/Nt*2*pi;
        t4 = t4/Nt*2*pi;
        t5 = t5/Nt*2*pi;
        t6 = t6/Nt*2*pi;
        t7 = t7/Nt*2*pi;
    end

    nd=1;
    zeta = [];
    tt = [];
    
    ncycle = 5;
    for k=1:ncycle

        if(ntransition(o)==2)
            if(k==1)
               a1=1;
               p1=0;
            else
               a1=a3;
               p1=p3;
            end
            
            clear p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7
            syms  p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7

            eq1 = a1*cos(o2*((k-1)*2*pi+t1)+p1)-a2*cos(o3*((k-1)*2*pi+t1)+p2)==0;
            eq2 = o2*a1*sin(o2*((k-1)*2*pi+t1)+p1)-o3*a2*sin(o3*((k-1)*2*pi+t1)+p2)==0;
            S = solve(eq1,eq2);
            p2 = double(S.p2(nd));
            a2 =double(S.a2(nd));
            
            eq3 = a2*cos(o3*((k-1)*2*pi+t2)+p2)-a3*cos(o2*((k-1)*2*pi+t2)+p3)==0;
            eq4 = o3*a2*sin(o3*((k-1)*2*pi+t2)+p2)-o2*a3*sin(o2*((k-1)*2*pi+t2)+p3)==0;
            S = solve(eq3,eq4);
            p3 = double(S.p3(nd));
            a3 =double(S.a3(nd));

            z1 = a1*cos(o2*((k-1)*2*pi+tt1)+p1);
            z2 = a2*cos(o3*((k-1)*2*pi+tt2)+p2);
            z3 = a3*cos(o2*((k-1)*2*pi+tt3)+p3);
            z4=[];
            z5=[];
            z6=[];
            z7=[];

            zeta = [zeta z1 z2 z3 z4 z5 z6 z7];
            
            tt = [tt (k-1)*2*pi+tt_hat];
            plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',blue)
            plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',orange)
            plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',blue)
            
        elseif(ntransition(o)==4)

            if(k==1)
               a1=1;
               p1=0;
            else
               a1=a5;
               p1=p5;
            end
            
            clear p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7
            syms  p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7

            eq1 = a1*cos(o2*((k-1)*2*pi+t1)+p1)-a2*cos(o1*((k-1)*2*pi+t1)+p2)==0;
            eq2 = o2*a1*sin(o2*((k-1)*2*pi+t1)+p1)-o1*a2*sin(o1*((k-1)*2*pi+t1)+p2)==0;
            S = solve(eq1,eq2);
            p2 = double(S.p2(nd));
            a2 =double(S.a2(nd));
            
            eq3 = a2*cos(o1*((k-1)*2*pi+t2)+p2)-a3*cos(o2*((k-1)*2*pi+t2)+p3)==0;
            eq4 = o1*a2*sin(o1*((k-1)*2*pi+t2)+p2)-o2*a3*sin(o2*((k-1)*2*pi+t2)+p3)==0;
            S = solve(eq3,eq4);
            p3 = double(S.p3(nd));
            a3 =double(S.a3(nd));
            
            eq5 = a3*cos(o2*((k-1)*2*pi+t3)+p3)-a4*cos(o3*((k-1)*2*pi+t3)+p4)==0;
            eq6 = o2*a3*sin(o2*((k-1)*2*pi+t3)+p3)-o3*a4*sin(o3*((k-1)*2*pi+t3)+p4)==0;
            S = solve(eq5,eq6);
            p4 = double(S.p4(nd));
            a4 =double(S.a4(nd));
            
            eq7 = a4*cos(o3*((k-1)*2*pi+t4)+p4)-a5*cos(o2*((k-1)*2*pi+t4)+p5)==0;
            eq8 = o3*a4*sin(o3*((k-1)*2*pi+t4)+p4)-o2*a5*sin(o2*((k-1)*2*pi+t4)+p5)==0;
            S = solve(eq7,eq8);
            p5 = double(S.p5(nd));
            a5 =double(S.a5(nd));

            z1 = a1*cos(o2*((k-1)*2*pi+tt1)+p1);
            z2 = a2*cos(o1*((k-1)*2*pi+tt2)+p2);
            z3 = a3*cos(o2*((k-1)*2*pi+tt3)+p3);
            z4 = a4*cos(o3*((k-1)*2*pi+tt4)+p4);
            z5 = a5*cos(o2*((k-1)*2*pi+tt5)+p5);
            z6=[];
            z7=[];

            zeta = [zeta z1 z2 z3 z4 z5 z6 z7];
            
            tt = [tt (k-1)*2*pi+tt_hat];
            plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',blue)
            plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',green)
            plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',blue)
            plot(k-1+tt4/2/pi,z4,'LineWidth',2,'Color',orange)
            plot(k-1+tt5/2/pi,z5,'LineWidth',2,'Color',blue)

        elseif(ntransition(o)==6)
            if(k==1)
               a1=1;
               p1=0;
            else
               a1=a7;
               p1=p7;
            end
            
            clear p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7
            syms  p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7

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
            
            eq7 = a4*cos(o2*((k-1)*2*pi+t4)+p4)-a5*cos(o3*((k-1)*2*pi+t4)+p5)==0;
            eq8 = o2*a4*sin(o2*((k-1)*2*pi+t4)+p4)-o3*a5*sin(o3*((k-1)*2*pi+t4)+p5)==0;
            S = solve(eq7,eq8);
            p5 = double(S.p5(nd));
            a5 =double(S.a5(nd));

            eq9 = a5*cos(o3*((k-1)*2*pi+t5)+p5)-a6*cos(o2*((k-1)*2*pi+t5)+p6)==0;
            eq10 = o3*a5*sin(o3*((k-1)*2*pi+t5)+p5)-o2*a6*sin(o2*((k-1)*2*pi+t5)+p6)==0;
            S = solve(eq9,eq10);
            p6 = double(S.p6(nd));
            a6 =double(S.a6(nd));
            
            eq11 = a6*cos(o2*((k-1)*2*pi+t6)+p6)-a7*cos(o1*((k-1)*2*pi+t6)+p7)==0;
            eq12 = o2*a6*sin(o2*((k-1)*2*pi+t6)+p6)-o1*a7*sin(o1*((k-1)*2*pi+t6)+p7)==0;
            S = solve(eq11,eq12);
            p7 = double(S.p7(nd));
            a7 =double(S.a7(nd));

            z1 = a1*cos(o1*((k-1)*2*pi+tt1)+p1);
            z2 = a2*cos(o2*((k-1)*2*pi+tt2)+p2);
            z3 = a3*cos(o1*((k-1)*2*pi+tt3)+p3);
            z4 = a4*cos(o2*((k-1)*2*pi+tt4)+p4);
            z5 = a5*cos(o3*((k-1)*2*pi+tt5)+p5);
            z6 = a6*cos(o2*((k-1)*2*pi+tt6)+p6);
            z7 = a7*cos(o1*((k-1)*2*pi+tt7)+p7);

            zeta = [zeta z1 z2 z3 z4 z5 z6 z7];
            
            tt = [tt (k-1)*2*pi+tt_hat];
            plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',green)
            plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',blue)
            plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',green)
            plot(k-1+tt4/2/pi,z4,'LineWidth',2,'Color',blue)
            plot(k-1+tt5/2/pi,z5,'LineWidth',2,'Color',orange)
            plot(k-1+tt6/2/pi,z6,'LineWidth',2,'Color',blue)
            plot(k-1+tt7/2/pi,z7,'LineWidth',2,'Color',green)

        end
        

    end

grid on;grid minor;set(gca,'FontSize',fontsize)
xt = 3.9;
yt = 18000;
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
xt = 3.9;
yt = 11;
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);

ylim([-14 13])



%%

xt = 0.78;
yt = 4;
str = '$R=2.4$';
o=97;

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


nexttile
plot(tt_hat/2/pi,ones(1,length(tt_hat)),'-','LineWidth',1,'Color',black);
hold on;
l1 = plot(tt_hat/2/pi,sigma1_imag,'LineWidth',2,'Color',brown);
l2 = plot(tt_hat/2/pi,sigma_new(o,:),'LineWidth',2,'Color',pink);
l3 = plot(tt_hat/2/pi,sigma1_imag_harmonic(o)*ones(1,length(tt_hat)),'--','LineWidth',2,'Color',yellow);
set(gca,'YScale', 'log');
ylim([1e-4 12])
grid on;grid minor;set(gca,'FontSize',fontsize)
title('Time-dependent frequency $\sigma(\hat t)$','Interpreter','latex','FontSize',fontsize+5);
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);






%%
nexttile
hold on;


    R = R_all(o);

    o1=om1(o);
    o2=om2(o);
    o3=om3(o);


    ntransition(o);
    transition_idx = find(diff(sigma_new_idx(o,:))~=0);

    if(ntransition(o)==2)
        t1 = transition_idx(1);
        t2 = transition_idx(2);
        t3 = Nt;
        tt1 = tt_hat(1:t1);
        tt2 = tt_hat(t1+1:t2);
        tt3 = tt_hat(t2+1:t3);
        t1 = t1/Nt*2*pi;
        t2 = t2/Nt*2*pi;
        t3 = t3/Nt*2*pi;
    elseif(ntransition(o)==4)
        t1 = transition_idx(1);
        t2 = transition_idx(2);
        t3 = transition_idx(3);
        t4 = transition_idx(4);
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
    elseif(ntransition(o)==6)
        t1 = transition_idx(1);
        t2 = transition_idx(2);
        t3 = transition_idx(3);
        t4 = transition_idx(4);
        t5 = transition_idx(5);
        t6 = transition_idx(6);
        t7 = Nt;
        tt1 = tt_hat(1:t1);
        tt2 = tt_hat(t1+1:t2);
        tt3 = tt_hat(t2+1:t3);
        tt4 = tt_hat(t3+1:t4);
        tt5 = tt_hat(t4+1:t5);
        tt6 = tt_hat(t5+1:t6);
        tt7 = tt_hat(t6+1:t7);
        t1 = t1/Nt*2*pi;
        t2 = t2/Nt*2*pi;
        t3 = t3/Nt*2*pi;
        t4 = t4/Nt*2*pi;
        t5 = t5/Nt*2*pi;
        t6 = t6/Nt*2*pi;
        t7 = t7/Nt*2*pi;
    end

    nd=1;
    zeta = [];
    tt = [];
    
    ncycle = 5;
    for k=1:ncycle

        if(ntransition(o)==2)
            if(k==1)
               a1=1;
               p1=0;
            else
               a1=a3;
               p1=p3;
            end
            
            clear p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7
            syms  p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7

            eq1 = a1*cos(o2*((k-1)*2*pi+t1)+p1)-a2*cos(o3*((k-1)*2*pi+t1)+p2)==0;
            eq2 = o2*a1*sin(o2*((k-1)*2*pi+t1)+p1)-o3*a2*sin(o3*((k-1)*2*pi+t1)+p2)==0;
            S = solve(eq1,eq2);
            p2 = double(S.p2(nd));
            a2 =double(S.a2(nd));
            
            eq3 = a2*cos(o3*((k-1)*2*pi+t2)+p2)-a3*cos(o2*((k-1)*2*pi+t2)+p3)==0;
            eq4 = o3*a2*sin(o3*((k-1)*2*pi+t2)+p2)-o2*a3*sin(o2*((k-1)*2*pi+t2)+p3)==0;
            S = solve(eq3,eq4);
            p3 = double(S.p3(nd));
            a3 =double(S.a3(nd));

            z1 = a1*cos(o2*((k-1)*2*pi+tt1)+p1);
            z2 = a2*cos(o3*((k-1)*2*pi+tt2)+p2);
            z3 = a3*cos(o2*((k-1)*2*pi+tt3)+p3);
            z4=[];
            z5=[];
            z6=[];
            z7=[];

            zeta = [zeta z1 z2 z3 z4 z5 z6 z7];
            
            tt = [tt (k-1)*2*pi+tt_hat];
            plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',blue)
            plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',orange)
            plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',blue)
            
        elseif(ntransition(o)==4)

            if(k==1)
               a1=1;
               p1=0;
            else
               a1=a5;
               p1=p5;
            end
            
            clear p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7
            syms  p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7

            eq1 = a1*cos(o2*((k-1)*2*pi+t1)+p1)-a2*cos(o1*((k-1)*2*pi+t1)+p2)==0;
            eq2 = o2*a1*sin(o2*((k-1)*2*pi+t1)+p1)-o1*a2*sin(o1*((k-1)*2*pi+t1)+p2)==0;
            S = solve(eq1,eq2);
            p2 = double(S.p2(nd));
            a2 =double(S.a2(nd));
            
            eq3 = a2*cos(o1*((k-1)*2*pi+t2)+p2)-a3*cos(o2*((k-1)*2*pi+t2)+p3)==0;
            eq4 = o1*a2*sin(o1*((k-1)*2*pi+t2)+p2)-o2*a3*sin(o2*((k-1)*2*pi+t2)+p3)==0;
            S = solve(eq3,eq4);
            p3 = double(S.p3(nd));
            a3 =double(S.a3(nd));
            
            eq5 = a3*cos(o2*((k-1)*2*pi+t3)+p3)-a4*cos(o3*((k-1)*2*pi+t3)+p4)==0;
            eq6 = o2*a3*sin(o2*((k-1)*2*pi+t3)+p3)-o3*a4*sin(o3*((k-1)*2*pi+t3)+p4)==0;
            S = solve(eq5,eq6);
            p4 = double(S.p4(nd));
            a4 =double(S.a4(nd));
            
            eq7 = a4*cos(o3*((k-1)*2*pi+t4)+p4)-a5*cos(o2*((k-1)*2*pi+t4)+p5)==0;
            eq8 = o3*a4*sin(o3*((k-1)*2*pi+t4)+p4)-o2*a5*sin(o2*((k-1)*2*pi+t4)+p5)==0;
            S = solve(eq7,eq8);
            p5 = double(S.p5(nd));
            a5 =double(S.a5(nd));

            z1 = a1*cos(o2*((k-1)*2*pi+tt1)+p1);
            z2 = a2*cos(o1*((k-1)*2*pi+tt2)+p2);
            z3 = a3*cos(o2*((k-1)*2*pi+tt3)+p3);
            z4 = a4*cos(o3*((k-1)*2*pi+tt4)+p4);
            z5 = a5*cos(o2*((k-1)*2*pi+tt5)+p5);
            z6=[];
            z7=[];

            zeta = [zeta z1 z2 z3 z4 z5 z6 z7];
            
            tt = [tt (k-1)*2*pi+tt_hat];
            plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',green)
            plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',blue)
            plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',orange)
            plot(k-1+tt4/2/pi,z4,'LineWidth',2,'Color',blue)
            plot(k-1+tt5/2/pi,z5,'LineWidth',2,'Color',green)

        elseif(ntransition(o)==6)
            if(k==1)
               a1=1;
               p1=0;
            else
               a1=a7;
               p1=p7;
            end
            
            clear p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7
            syms  p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7

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
            
            eq7 = a4*cos(o2*((k-1)*2*pi+t4)+p4)-a5*cos(o3*((k-1)*2*pi+t4)+p5)==0;
            eq8 = o2*a4*sin(o2*((k-1)*2*pi+t4)+p4)-o3*a5*sin(o3*((k-1)*2*pi+t4)+p5)==0;
            S = solve(eq7,eq8);
            p5 = double(S.p5(nd));
            a5 =double(S.a5(nd));

            eq9 = a5*cos(o3*((k-1)*2*pi+t5)+p5)-a6*cos(o2*((k-1)*2*pi+t5)+p6)==0;
            eq10 = o3*a5*sin(o3*((k-1)*2*pi+t5)+p5)-o2*a6*sin(o2*((k-1)*2*pi+t5)+p6)==0;
            S = solve(eq9,eq10);
            p6 = double(S.p6(nd));
            a6 =double(S.a6(nd));
            
            eq11 = a6*cos(o2*((k-1)*2*pi+t6)+p6)-a7*cos(o1*((k-1)*2*pi+t6)+p7)==0;
            eq12 = o2*a6*sin(o2*((k-1)*2*pi+t6)+p6)-o1*a7*sin(o1*((k-1)*2*pi+t6)+p7)==0;
            S = solve(eq11,eq12);
            p7 = double(S.p7(nd));
            a7 =double(S.a7(nd));

            z1 = a1*cos(o1*((k-1)*2*pi+tt1)+p1);
            z2 = a2*cos(o2*((k-1)*2*pi+tt2)+p2);
            z3 = a3*cos(o1*((k-1)*2*pi+tt3)+p3);
            z4 = a4*cos(o2*((k-1)*2*pi+tt4)+p4);
            z5 = a5*cos(o3*((k-1)*2*pi+tt5)+p5);
            z6 = a6*cos(o2*((k-1)*2*pi+tt6)+p6);
            z7 = a7*cos(o1*((k-1)*2*pi+tt7)+p7);

            zeta = [zeta z1 z2 z3 z4 z5 z6 z7];
            
            tt = [tt (k-1)*2*pi+tt_hat];
            plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',green)
            plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',blue)
            plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',green)
            plot(k-1+tt4/2/pi,z4,'LineWidth',2,'Color',blue)
            plot(k-1+tt5/2/pi,z5,'LineWidth',2,'Color',orange)
            plot(k-1+tt6/2/pi,z6,'LineWidth',2,'Color',blue)
            plot(k-1+tt7/2/pi,z7,'LineWidth',2,'Color',green)

        end
        

    end



grid on;grid minor;set(gca,'FontSize',fontsize)
xt = 3.9;
yt = 1500;
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
xt = 3.9;
yt = 11;
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);
ylim([-14 13])





%%

xt = 0.78;
yt = 4;
str = '$R=3.3$';
o=133;

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


nexttile
plot(tt_hat/2/pi,ones(1,length(tt_hat)),'-','LineWidth',1,'Color',black);
hold on;
l1 = plot(tt_hat/2/pi,sigma1_imag,'LineWidth',2,'Color',brown);
l2 = plot(tt_hat/2/pi,sigma_new(o,:),'LineWidth',2,'Color',pink);
l3 = plot(tt_hat/2/pi,sigma1_imag_harmonic(o)*ones(1,length(tt_hat)),'--','LineWidth',2,'Color',yellow);
set(gca,'YScale', 'log');
ylim([1e-4 12])
grid on;grid minor;set(gca,'FontSize',fontsize)
title('Time-dependent frequency $\sigma(\hat t)$','Interpreter','latex','FontSize',fontsize+5);
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);


xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');









%%
nexttile
hold on;


    R = R_all(o);

    o1=om1(o);
    o2=om2(o);
    o3=om3(o);


    ntransition(o);
    transition_idx = find(diff(sigma_new_idx(o,:))~=0);

    if(ntransition(o)==2)
        t1 = transition_idx(1);
        t2 = transition_idx(2);
        t3 = Nt;
        tt1 = tt_hat(1:t1);
        tt2 = tt_hat(t1+1:t2);
        tt3 = tt_hat(t2+1:t3);
        t1 = t1/Nt*2*pi;
        t2 = t2/Nt*2*pi;
        t3 = t3/Nt*2*pi;
    elseif(ntransition(o)==4)
        t1 = transition_idx(1);
        t2 = transition_idx(2);
        t3 = transition_idx(3);
        t4 = transition_idx(4);
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
    elseif(ntransition(o)==6)
        t1 = transition_idx(1);
        t2 = transition_idx(2);
        t3 = transition_idx(3);
        t4 = transition_idx(4);
        t5 = transition_idx(5);
        t6 = transition_idx(6);
        t7 = Nt;
        tt1 = tt_hat(1:t1);
        tt2 = tt_hat(t1+1:t2);
        tt3 = tt_hat(t2+1:t3);
        tt4 = tt_hat(t3+1:t4);
        tt5 = tt_hat(t4+1:t5);
        tt6 = tt_hat(t5+1:t6);
        tt7 = tt_hat(t6+1:t7);
        t1 = t1/Nt*2*pi;
        t2 = t2/Nt*2*pi;
        t3 = t3/Nt*2*pi;
        t4 = t4/Nt*2*pi;
        t5 = t5/Nt*2*pi;
        t6 = t6/Nt*2*pi;
        t7 = t7/Nt*2*pi;
    end

    nd=1;
    zeta = [];
    tt = [];
    
    ncycle = 5;
    for k=1:ncycle

        if(ntransition(o)==2)
            if(k==1)
               a1=1;
               p1=0;
            else
               a1=a3;
               p1=p3;
            end
            
            clear p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7
            syms  p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7

            eq1 = a1*cos(o2*((k-1)*2*pi+t1)+p1)-a2*cos(o3*((k-1)*2*pi+t1)+p2)==0;
            eq2 = o2*a1*sin(o2*((k-1)*2*pi+t1)+p1)-o3*a2*sin(o3*((k-1)*2*pi+t1)+p2)==0;
            S = solve(eq1,eq2);
            p2 = double(S.p2(nd));
            a2 =double(S.a2(nd));
            
            eq3 = a2*cos(o3*((k-1)*2*pi+t2)+p2)-a3*cos(o2*((k-1)*2*pi+t2)+p3)==0;
            eq4 = o3*a2*sin(o3*((k-1)*2*pi+t2)+p2)-o2*a3*sin(o2*((k-1)*2*pi+t2)+p3)==0;
            S = solve(eq3,eq4);
            p3 = double(S.p3(nd));
            a3 =double(S.a3(nd));

            z1 = a1*cos(o2*((k-1)*2*pi+tt1)+p1);
            z2 = a2*cos(o3*((k-1)*2*pi+tt2)+p2);
            z3 = a3*cos(o2*((k-1)*2*pi+tt3)+p3);
            z4=[];
            z5=[];
            z6=[];
            z7=[];

            zeta = [zeta z1 z2 z3 z4 z5 z6 z7];
            
            tt = [tt (k-1)*2*pi+tt_hat];
            plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',blue)
            plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',orange)
            plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',blue)
            
        elseif(ntransition(o)==4)

            if(k==1)
               a1=1;
               p1=0;
            else
               a1=a5;
               p1=p5;
            end
            
            clear p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7
            syms  p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7

            eq1 = a1*cos(o2*((k-1)*2*pi+t1)+p1)-a2*cos(o1*((k-1)*2*pi+t1)+p2)==0;
            eq2 = o2*a1*sin(o2*((k-1)*2*pi+t1)+p1)-o1*a2*sin(o1*((k-1)*2*pi+t1)+p2)==0;
            S = solve(eq1,eq2);
            p2 = double(S.p2(nd));
            a2 =double(S.a2(nd));
            
            eq3 = a2*cos(o1*((k-1)*2*pi+t2)+p2)-a3*cos(o2*((k-1)*2*pi+t2)+p3)==0;
            eq4 = o1*a2*sin(o1*((k-1)*2*pi+t2)+p2)-o2*a3*sin(o2*((k-1)*2*pi+t2)+p3)==0;
            S = solve(eq3,eq4);
            p3 = double(S.p3(nd));
            a3 =double(S.a3(nd));
            
            eq5 = a3*cos(o2*((k-1)*2*pi+t3)+p3)-a4*cos(o3*((k-1)*2*pi+t3)+p4)==0;
            eq6 = o2*a3*sin(o2*((k-1)*2*pi+t3)+p3)-o3*a4*sin(o3*((k-1)*2*pi+t3)+p4)==0;
            S = solve(eq5,eq6);
            p4 = double(S.p4(nd));
            a4 =double(S.a4(nd));
            
            eq7 = a4*cos(o3*((k-1)*2*pi+t4)+p4)-a5*cos(o2*((k-1)*2*pi+t4)+p5)==0;
            eq8 = o3*a4*sin(o3*((k-1)*2*pi+t4)+p4)-o2*a5*sin(o2*((k-1)*2*pi+t4)+p5)==0;
            S = solve(eq7,eq8);
            p5 = double(S.p5(nd));
            a5 =double(S.a5(nd));

            z1 = a1*cos(o2*((k-1)*2*pi+tt1)+p1);
            z2 = a2*cos(o1*((k-1)*2*pi+tt2)+p2);
            z3 = a3*cos(o2*((k-1)*2*pi+tt3)+p3);
            z4 = a4*cos(o3*((k-1)*2*pi+tt4)+p4);
            z5 = a5*cos(o2*((k-1)*2*pi+tt5)+p5);
            z6=[];
            z7=[];

            zeta = [zeta z1 z2 z3 z4 z5 z6 z7];
            
            tt = [tt (k-1)*2*pi+tt_hat];
            plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',green)
            plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',blue)
            plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',orange)
            plot(k-1+tt4/2/pi,z4,'LineWidth',2,'Color',blue)
            plot(k-1+tt5/2/pi,z5,'LineWidth',2,'Color',green)

        elseif(ntransition(o)==6)
            if(k==1)
               a1=1;
               p1=0;
            else
               a1=a7;
               p1=p7;
            end
            
            clear p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7
            syms  p2 a2 p3 a3 p4 a4 p5 a5 a6 p6 a7 p7

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
            
            eq7 = a4*cos(o2*((k-1)*2*pi+t4)+p4)-a5*cos(o3*((k-1)*2*pi+t4)+p5)==0;
            eq8 = o2*a4*sin(o2*((k-1)*2*pi+t4)+p4)-o3*a5*sin(o3*((k-1)*2*pi+t4)+p5)==0;
            S = solve(eq7,eq8);
            p5 = double(S.p5(nd));
            a5 =double(S.a5(nd));

            eq9 = a5*cos(o3*((k-1)*2*pi+t5)+p5)-a6*cos(o2*((k-1)*2*pi+t5)+p6)==0;
            eq10 = o3*a5*sin(o3*((k-1)*2*pi+t5)+p5)-o2*a6*sin(o2*((k-1)*2*pi+t5)+p6)==0;
            S = solve(eq9,eq10);
            p6 = double(S.p6(nd));
            a6 =double(S.a6(nd));
            
            eq11 = a6*cos(o2*((k-1)*2*pi+t6)+p6)-a7*cos(o1*((k-1)*2*pi+t6)+p7)==0;
            eq12 = o2*a6*sin(o2*((k-1)*2*pi+t6)+p6)-o1*a7*sin(o1*((k-1)*2*pi+t6)+p7)==0;
            S = solve(eq11,eq12);
            p7 = double(S.p7(nd));
            a7 =double(S.a7(nd));

            z1 = a1*cos(o1*((k-1)*2*pi+tt1)+p1);
            z2 = a2*cos(o2*((k-1)*2*pi+tt2)+p2);
            z3 = a3*cos(o1*((k-1)*2*pi+tt3)+p3);
            z4 = a4*cos(o2*((k-1)*2*pi+tt4)+p4);
            z5 = a5*cos(o3*((k-1)*2*pi+tt5)+p5);
            z6 = a6*cos(o2*((k-1)*2*pi+tt6)+p6);
            z7 = a7*cos(o1*((k-1)*2*pi+tt7)+p7);

            zeta = [zeta z1 z2 z3 z4 z5 z6 z7];
            
            tt = [tt (k-1)*2*pi+tt_hat];
            plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',green)
            plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',blue)
            plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',green)
            plot(k-1+tt4/2/pi,z4,'LineWidth',2,'Color',blue)
            plot(k-1+tt5/2/pi,z5,'LineWidth',2,'Color',orange)
            plot(k-1+tt6/2/pi,z6,'LineWidth',2,'Color',blue)
            plot(k-1+tt7/2/pi,z7,'LineWidth',2,'Color',green)

        end
        

    end



grid on;grid minor;set(gca,'FontSize',fontsize)
xt = 3.9;
yt = 17.5e4;
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
yt = 11;
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+4);
xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');
ylim([-14 13])





tiledlay.TileSpacing = 'compact';
tiledlay.Padding = 'compact';




AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal')



print('-dpng','-r300',['fig_supp_new/figS_ana_topo1_matlab.png']);
