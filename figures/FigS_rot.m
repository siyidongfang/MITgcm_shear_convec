clear;
addpath ../instability_rotating/
addpath ../instability_rotating/exps_rotating/

fname = 'flat_N1e-3output';
% fname = 'noV_flat_N1e-3output';
strn = 'Flat bottom';


load([fname '.mat'])
addpath ../analysis/colormaps/
% grow_rw = grow_rw_ke;
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
set(gcf,'Color','w','Position', [-46 153 1339 824])
subplot(2,2,1)
pcolor(R,1./rw_all,grow_rw');
shading flat;colorbar;
colormap(WhiteBlueGreenYellowRed(0))
% colormap([[1 1 1]*0.97;jet])
set(gca,'FontSize',fontsize)
xlabel('${R_i}_\mathrm{min}^{-1}$','Interpreter','latex')
ylabel('$m_0/k_0$','Interpreter','latex')
title(strn,'Interpreter','latex','FontSize',fontsize+5)
% clim([0 0.4])
% xlim([0 10])
% ylim([-25 25])
set(gca,'TickDir','out');
grid on;grid minor;





subplot(2,2,2)
plot(R,max(grow_rw,[],2),'LineWidth',2);
hold on;
% load('../instability_km/exps_new/topo4_nu0_output.mat','rw_all','grow_rw','shear_all')
% load('../figures/fig4/Ri_topo4.mat')
load('../instability_km/exps_varyingN/N1e-3output','grow_rw','shear_all')
load('../figures/fig4/Ri_flat.mat')
max_grow_km = max(grow_rw,[],2);
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end
plot(1./Ri_km,max_grow_km,'LineWidth',2);
hold off;
grid on;grid minor;
set(gca,'FontSize',fontsize)
axis tight;
xlim([0 6.25])
xlabel('${R_i}_\mathrm{min}^{-1}$','Interpreter','latex','FontSize',fontsize+4)
ylabel('Maximum growth rate (1/hour)','Interpreter','latex')
title(strn,'Interpreter','latex','FontSize',fontsize+5)
legend('Rotating','Non-rotating','Position',[0.6004 0.8076 0.0993 0.0516],'FontSize',fontsize+3,'Interpreter','latex');
legend('Boxoff')

set(gca,'TickDir','out');




fname = 'topo4_N1e-3output';
% fname = 'noV_topo4_N1e-3output';
strn = 'Sloping bottom';


load([fname '.mat'])
addpath ../analysis/colormaps/
% grow_rw = grow_rw_ke;
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


set(gcf,'Color','w','Position', [-46 153 1339 824])
subplot(2,2,3)
pcolor(R,1./rw_all,grow_rw');
shading flat;colorbar;
% colormap(WhiteBlueGreenYellowRed(0))
set(gca,'FontSize',fontsize)
xlabel('${R_i}_\mathrm{min}^{-1}$','Interpreter','latex')
ylabel('$m_0/k_0$','Interpreter','latex')
title(strn,'Interpreter','latex','FontSize',fontsize+5)
% clim([0 0.4])
% xlim([0 10])
% ylim([-25 25])
set(gca,'TickDir','out');
grid on;grid minor;





subplot(2,2,4)
plot(R,max(grow_rw,[],2),'LineWidth',2);
hold on;
% load('../instability_km/exps_new/topo4_nu0_output.mat','rw_all','grow_rw','shear_all')
% load('../figures/fig4/Ri_topo4.mat')
load('../instability_km/exps_varyingN/N1e-3output','grow_rw','shear_all')
load('../figures/fig4/Ri_flat.mat')
max_grow_km = max(grow_rw,[],2);
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end
plot(1./Ri_km,max_grow_km,'LineWidth',2);
hold off;
grid on;grid minor;
set(gca,'FontSize',fontsize)
axis tight;
xlim([0 5.1])
xlabel('${R_i}_\mathrm{min}^{-1}$','Interpreter','latex','FontSize',fontsize+4)
ylabel('Maximum growth rate (1/hour)','Interpreter','latex')
title(strn,'Interpreter','latex','FontSize',fontsize+5)
l2 = legend('Rotating','Non-rotating','Position',[0.6019 0.3331 0.0993 0.0516],'FontSize',fontsize+3,'Interpreter','latex');
legend('Boxoff')
set(gca,'TickDir','out');




AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal')




%%
print('-dpng','-r300',['fig_supp_new/figS_rot2_matlab.png']);

