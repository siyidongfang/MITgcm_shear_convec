clear;

addpath ../analysis/functions
load_colors

load('freq_topo_3om.mat')

% for o = 1:177 
for o = 2:241 
    o
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

    
    % figure(1)
    % clf;set(gcf,'Color','w','Position',[41 146 1275 428])
    % hold on;
    
    nd=1;
    zeta = [];
    tt = [];
    
    ncycle = 40;
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
            % plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',blue)
            % plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',orange)
            % plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',blue)
            
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
            % plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',green)
            % plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',blue)
            % plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',orange)
            % plot(k-1+tt4/2/pi,z4,'LineWidth',2,'Color',blue)
            % plot(k-1+tt5/2/pi,z5,'LineWidth',2,'Color',green)

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
            % plot(k-1+tt1/2/pi,z1,'LineWidth',2,'Color',green)
            % plot(k-1+tt2/2/pi,z2,'LineWidth',2,'Color',blue)
            % plot(k-1+tt3/2/pi,z3,'LineWidth',2,'Color',green)
            % plot(k-1+tt4/2/pi,z4,'LineWidth',2,'Color',blue)
            % plot(k-1+tt5/2/pi,z5,'LineWidth',2,'Color',orange)
            % plot(k-1+tt6/2/pi,z6,'LineWidth',2,'Color',blue)
            % plot(k-1+tt7/2/pi,z7,'LineWidth',2,'Color',green)

        end
        

    end


% % ylim([-0.6 1.6]*1e4)
% grid on;grid minor;set(gca,'Fontsize',20);
% xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');
% title(['$\zeta(\hat t)$, R=' num2str(R,3)],'Interpreter','latex');
% print('-dpng','-r150',['figures_topo4/topo4_o' num2str(o) '.png']);

ke = 0.5*zeta.^2;
fit_span = 1:Nt*ncycle;
xxplot = tt/2/pi*12; %%% in hours
yyplot = log(ke)/2;
[pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
grow(o)=pp(1);
[y_fit,delta_fit] = polyval(pp,xxplot,S);

% figure(2)
% clf;set(gcf,'Color','w','Position',[41 146 1275 428])
% plot(tt/2/pi,log(ke)/2,'LineWidth',2)
% hold on;
% plot(xxplot(fit_span)/12, y_fit(fit_span));
% hold off;
% grid on;grid minor;set(gca,'Fontsize',20);
% ylabel('$\ln(\mathrm{TKE})/2$','Interpreter','latex');
% xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');
% title(['TKE$\sim0.5*\zeta^2(\hat t)$, R=' num2str(R,3)],'Interpreter','latex');
% print('-dpng','-r150',['figures_topo4/ke' num2str(o) '.png']);


end

close all;
save('grow_topo_3om.mat')


%%
load('grow_topo_3om.mat')
figure(3)
clf;set(gcf,'Color','w','Position',[41 146 700 428])
plot(R_all,grow,'LineWidth',2)
grid on;grid minor;set(gca,'Fontsize',20);
title('Growth rate (1/hour)','Interpreter','latex');


% load('/Users/ysi/MITgcm_shear_convec/instability_km/exps_new/topo4_nu0_output.mat','grow_rw','shear_all')
load('../instability_km/exps_new/topo4_nu0_output.mat','rw_all','grow_rw','shear_all')
load('../figures/fig4/Ri_topo4.mat')
hold on;

max_grow_km = max(grow_rw,[],2);
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));

    r=cosd(topo)/sind(topo)-N_hat*sqrt(1/Ri_km(i))

    [a(i) b(i)] = min(abs(r-1./rw_all));
    grow_mzero(i) = grow_rw(i,b(i));
end

plot(1./Ri_km,grow_mzero,'--','LineWidth',2)
xlim([0 4.4])
ylim([-0.04 0.35])
xlabel('R = ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
legend('Analytical solution','Theory')


print('-dpng','-r150',['figures_topo4/grow_topo.png']);
