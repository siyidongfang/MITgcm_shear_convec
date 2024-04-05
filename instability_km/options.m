%%%%% All variables are dimensional variables

clear; close all;

NOdiff = true;
nt_percycle = 72; %%% nt_percycle = 720 if use diffusion.
NTtide = 100;
topo = 0;
N = sqrt(10)*1e-3;
omega = 0.1*1e-3;
Ptide = 2*pi/omega;
m0 =1;
shear = sqrt(N^2);
rs = shear/omega; %%% shear over omega 
constants;

rw_all = 0.001:0.0005:0.01;

for i=1:length(rw_all)
    rw = rw_all(i);
    kx = m0*rw;
    loop;
    if(grow(i)>0 && grow(i)<1e-3)
        NTtide = 300;
        constants;
        loop;
    end
end





% N = 1e-3;
% shear = 9.7e-4;
% shear_all = [0:0.1:3]*1e-3;
% ns = length(shear_all);
% Ptide = 43200;
% omega = 2*pi/Ptide;
% topo_all = [0:9];
% mz_all = [0:0.09:6];
% kx_all = [-0.1:0.003:0.1];
% rw_all = 10.^([-2:0.1:-1.2 -1.15:0.05:0.6 1 2 3 4]); %%% kx/mz
% rw_all = 10.^([-3:0.005:1]); %%% kx/mz
% rw_all = 10.^([-3:0.02:1]); %%% kx/mz
% nr = length(rw_all);


        % growth(k,m) =  pKE(1);
        % save(['output_Ri1_Nsq1e-6_topo4/mz' num2str(mz) 'kx' num2str(kx) 'shear' num2str(shear) '.mat'])
    % end
% end

% save(['output_topo0/growth_topo' num2str(topo) '_mz' num2str(mz) 'kx' num2str(kx) '_NOdiff.mat'])

% end

% if(nr>1)
% figure(21);
% % clf
% set(gcf,'Color','w');
% plot(log10(rw_all),growth,'LineWidth',2)
% hold on;
% grid on;grid minor;
% title('Growth rate (1/hour)')
% ylabel('(1/hour)')
% xlabel('Wavenumber ratio log10(k_x/m_z)')
% set(gca,'fontsize',20)
% end





