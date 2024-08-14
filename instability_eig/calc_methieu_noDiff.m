%%
clear all;addpath ../analysis/colormaps/

expdir = 'shear_methieu_noDiff/'
Ptide = 43200;
topo = 0;
N = 1e-3;
Ns = 1001;
omega = 2*pi/Ptide;
shear_all = [0:3*N/(Ns-1):3*N]; 

grow_all = [];
omega0_all = [];
epsilon_all = [];
Ri_all = [];

h_shear=250;
lam_x_all = [-8000:1:8000];
m0_rw = 2*pi/h_shear;
kx_all = 2*pi./lam_x_all;
omega0 = abs(N*kx_all/m0_rw);

dt_ri = 10;
Nt_ri = 100*Ptide/dt_ri;
tt_ri = dt_ri:dt_ri:Nt_ri*dt_ri;

Ri_shear = NaN.*zeros(1,Ns); %%% minimum Ri as a function of shear

for ns=1:Ns
    shear = shear_all(ns);
    load([expdir 'ptide' num2str(Ptide) '_topo' num2str(topo) '_N' num2str(N*1e3) '_shear' num2str(shear*1e3) '.mat']);
    grow(grow<=0)=NaN;
    grow = log(grow)/43200*3600;
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

    Ri_all = [Ri_all Ri_shear(ns).*ones(1,Nk)];

end


[omega0_sort,I] = sort(omega0_all);
epsilon_sort = epsilon_all(I);
grow_sort = grow_all(I);
Ri_sort = Ri_all(I);

save('../figures/fig5/fig5_topo0_noDiff_eig_new.mat')


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
