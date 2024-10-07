%%
clear all;addpath ../analysis/colormaps/

load('products/grow_flat_eig_km.mat')
Nsidx = 1:759;

% load('products/grow_topo4_eig_km.mat')
% Nsidx = 1:756;


% figure(20)
% pcolor(grow_eig_km);shading flat;colorbar;colormap(WhiteBlueGreenYellowRed(0))
% clim([0 0.4])

omega0 = abs(N*kx_all(round(Nk/2):end)/m0_rw);

dt_ri = 10;
Nt_ri = 100*Ptide/dt_ri;
tt_ri = dt_ri:dt_ri:Nt_ri*dt_ri;

Ri_shear = NaN.*zeros(1,Ns); %%% minimum Ri as a function of shear
grow_all = [];
epsilon_all = [];
omega0_all = [];
Ri_all = [];

for ns=Nsidx
    shear = shear_all(ns);
    grow = grow_eig_km(round(Nk/2):end,ns);
    grow = grow';
    grow_all = [grow_all grow];

    epsilon = 2*omega0.^3/omega*shear/N;
    epsilon_all=[epsilon_all epsilon];
    omega0_all = [omega0_all omega0];

    %%% Compute Richardson number
    Ri_inverse = (shear*cos(omega*tt_ri)).^2./(N^2*cosd(topo) - N^2*sind(topo)/omega*shear*sin(omega*tt_ri));
    Ri_shear(ns) = 1/max(Ri_inverse);  
    if(min(Ri_inverse)<0)
        isConvec = 1
    end
    Ri_all = [Ri_all Ri_shear(ns).*ones(1,length(grow))];
end


[omega0_sort,I] = sort(omega0_all);
epsilon_sort = epsilon_all(I);
grow_sort = grow_all(I);
Ri_sort = Ri_all(I);

save('../figures/fig5/fig5_topo0_Ri_eig.mat')


%%

% grow_sort(grow_sort<0)=NaN;
% 
% figure(1)
% clf;set(gcf,'Color','w');
% scatter(omega0_sort/omega,epsilon_sort/omega^2,20,grow_sort,'filled')
% % scatter(1./Ri_sort,epsilon_sort/omega^2,20,grow_sort,'filled')
% colorbar;
% clim([0 0.3])
% % colormap(cmocean('balance'))
% colormap(jet)
% ylabel('$\sqrt\epsilon/\omega = \sqrt {2\Lambda\omega_0^3/(N\omega^3)}$','Interpreter','latex')
% title('Growth rate (1/hour)')
% xlabel('$\sqrt\delta/\omega = \omega_0 /\omega = N k_x/m_z/\omega$','Interpreter','latex')
% set(gca,'fontsize',20)
% xlim([0 2.5])
% % xlim([0 6])
% ylim([0 4])
% 
% 
% 
% 
% 
% 
% 
% 
% 
