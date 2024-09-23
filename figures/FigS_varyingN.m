clear;
addpath ../instability_km/exps_varyingN/

% fname = 'N1e-3output';
% strn = '$\tilde N=6.9\ \omega$';
% figname = 'figS_N0_matlab'
% 
% fname = 'N3omegaoutput';
% strn = '$\tilde N=3\ \omega$';
% figname = 'figS_N1_matlab'
% 
% fname = 'N20omegaoutput';
% strn = '$\tilde N=20\ \omega$';
% figname = 'figS_N2_matlab'
% 
fname = 'N100omegaoutput';
strn = '$\tilde N=100\ \omega$';
figname = 'figS_N3_matlab'

load([fname '.mat'])
addpath ../analysis/colormaps/
grow_rw(grow_rw<0)=NaN;

R = shear_all.^2/N^2;


grow_all = [];
omega0_all = [];
epsilon_all = [];
Ri_all = [];
omega0 = N*kx_all/m0_rw;

dt_ri = 10;
Nt_ri = 100*Ptide/dt_ri;
tt_ri = dt_ri:dt_ri:Nt_ri*dt_ri;

Ri_shear = NaN.*zeros(1,Ns); %%% minimum Ri as a function of shear

for ns=1:Ns
    shear = shear_all(ns);

    epsilon = 2*omega0.^3/omega*shear/N;
    epsilon_all=[epsilon_all epsilon];

    grow = grow_rw(ns,:);
    grow_all = [grow_all grow];

    omega0_all = [omega0_all omega0];


    %%% Compute Richardson number
    Ri_inverse = (shear*cos(omega*tt_ri)).^2./(N^2*cosd(topo) - N^2*sind(topo)/omega*shear*sin(omega*tt_ri));
    % figure(50);plot(tt_ri,Ri_inverse)
    Ri_shear(ns) = 1/max(Ri_inverse);  
    if(min(Ri_inverse)<0)
        isConvec = 1
    end

    Ri_all = [Ri_all Ri_shear(ns).*ones(1,Nk)];

end

[omega0_sort,I] = sort(omega0_all);
epsilon_sort = epsilon_all(I);
grow_sort = grow_all(I);
Ri_sort = Ri_all(I);



fontsize = 18;

fg1 = figure(1);
clf;
set(gcf,'Color','w','Position', [-46 153 1340 824])
subplot(2,2,1)
pcolor(R,1./rw_all,grow_rw');
shading flat;colorbar;
set(gca,'FontSize',fontsize)
xlabel('${R_i}_\mathrm{min}^{-1}$','Interpreter','latex')
ylabel('$m_0/k_0$','Interpreter','latex')
title(strn,'Interpreter','latex','FontSize',fontsize+5)
% clim([0 0.4])
xlim([0 9])
% ylim([0 450])
set(gca,'TickDir','out');
grid on;grid minor;



subplot(2,2,2)
plot(R,max(grow_rw,[],2),'LineWidth',2);
grid on;grid minor;
set(gca,'FontSize',fontsize)
axis tight;
xlim([0 9])
xlabel('${R_i}_\mathrm{min}^{-1}$','Interpreter','latex','FontSize',fontsize+4)
ylabel('Maximum growth rate (1/hour)','Interpreter','latex')
title(strn,'Interpreter','latex','FontSize',fontsize+5)
set(gca,'TickDir','out');

subplot(2,2,3)
scatter(omega0_sort/omega,epsilon_sort/omega^2,10,grow_sort,'filled')
shading flat;colorbar;
set(gca,'FontSize',fontsize)
title(strn,'Interpreter','latex','FontSize',fontsize+5)
xlim([0 2])
ylim([0 10])
% clim([0 max(grow_rw,[],'all')/3])
grid on;grid minor;
mycolormap=jet;
colormap([[1 1 1]*0.97;mycolormap])
% colormap(WhiteBlueGreenYellowRed(0))
set(gca,'TickDir','out');


ylabel('$\sqrt\epsilon_0/\omega = \sqrt {\frac{2\omega_0^3}{\omega}\frac{\Lambda}{\tilde N}}\Big/\omega$','Interpreter','latex','FontSize',fontsize+4)
xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);


subplot(2,2,4)
scatter(omega0_sort/omega,sqrt(epsilon_sort)/omega,20,grow_sort,'o','filled');
clim([0 0.4])
scatter(omega0_sort/omega,Ri_sort,10,grow_sort,'o', 'filled');
set(gca,'fontsize',fontsize)
ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
box on;
grid on;grid minor;colorbar;
xlim([0 2])
ylim([0.01 3000])
% clim([0 max(grow_rw,[],'all')/3])
set(gca, 'YScale', 'log')
set(gca,'TickDir','out');
title(strn,'Interpreter','latex','FontSize',fontsize+5)


% AddLetters2Plots(fg1, 'HShift', 0, 'VShift', 0, 'Direction', 'LeftRight')
AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal')

print('-dpng','-r300',['fig_supp_new/' figname '.png']);

