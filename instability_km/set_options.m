

% %%%%%% exps_topo4_diff %%%%%%
% expdir = 'exps_topo4_diff/'; 
% topo=4;
% N = sqrt(1)*1e-3;
% Ptide = 43200;
% omega = 2*pi/Ptide;
% shear_Ri0_25 = 0.0018;
% shear_Ri1 = 0.97e-3;
% shear_all = [0:1e-4:shear_Ri0_25];
% m0_all = [0:0.1:10];
% kx_all = [-0.5:0.0025*4:0.5];

% %%%%%% exps_KH %%%%%%
% expdir = 'exps_KH/';
% topo = 0;
% omega = 0;
% N = sqrt(1)*1e-3;
% Ptide = 43200; 
% shear_Ri0_25 = 2*N;
% shear_all = shear_Ri0_25*10;
% m0 =1;
% rw_all = 10.^([-2:0.1:5]); %%% For K-H instability

% %%%%%% exps_Radko %%%%%%
% expdir = 'exps_Radko/';
% topo = 0;
% N = sqrt(10)*1e-3;
% omega = 0.1*1e-3;
% Ptide = 2*pi/omega;
% shear_Ri0_25 = 0.0050; %%% calculated by the script calc_shear_from_Ri
% shear_Ri1 = 0.0031625;
% shear_Ri5 = 0.001415;
% shear_all = shear_Ri5;
% % m0_all = [0:0.015:6];
% % kx_all = [-0.1:0.0005:0.1];
% m0_all = [0:0.01:4];
% kx_all = [-0.5:0.0025:0.5];

% %%%%%% exps_topo4 %%%%%%
% expdir = 'exps_topo4/';
% topo=4;
% N = sqrt(1)*1e-3;
% Ptide = 43200;
% omega = 2*pi/Ptide;
% shear_Ri0_25 = 0.0018;
% shear_Ri1 = 0.97e-3;
% % load('rw_mg.mat') %%% The wavenumber ratio m0/k that corresponds to the largest growth rate with out diffusion; or load('rw_mg_test3.mat') with diffusion kappa=1e-5,nu=1e-6.
% % shear_all = [0:0.5e-5:shear_Ri0_25]; % Ri=1, shear = 0.97e-3;
% shear_all = [0:0.005:0.495]*1e-3;
% omega0_all = sqrt([0:0.01:6])*omega;
% rw_all = omega0_all/N;
% % rw_all= 10.^([-2:0.1:-1.1 -1.05 -1:0.0025:-0.7 -0.6:0.1:1]);
% m0 =1;
% load("rw_mg_exps_topo4.mat",'rw_mg')

% %%%%%% exps_flatbottom %%%%%%
% expdir = 'exps_flat/';
% topo=0;
% N = sqrt(1)*1e-3;
% Ptide = 43200;
% omega = 2*pi/Ptide;
% shear_Ri0_25 = 2*N;
% shear_Ri1 = N;
% omega0_all = sqrt([0:0.01:6])*omega;
% rw_all = omega0_all/N;
% % shear_all = [0:0.5e-5:shear_Ri0_25];
% shear_all = [0:0.005:0.355]*1e-3;
% 
% m0 =1;
% load("rw_mg_exps_flat.mat",'rw_mg')







% rw_all = 10.^([-1.5:0.1:-0.5]); 
% rw_all= 10.^([-3:0.1:-1 -0.95:0.01:0 0.1:0.1:1]);
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
% save(['output_topo0/growth_topo' num2str(topo) '_mz' num2str(mz) 'kx' num2str(kx) '_NOdiff.mat'])


% rw_all= 10.^([-2:0.1:-1 -0.95:0.01:-0.5 0.6:0.1:1]);
% rw_all= 10.^([-2:0.1:-1.1 -1:0.005:-0.7 -0.6:0.1:1]);
% rw_all = 10.^([-1.2:0.01:-0.5]);
% rw_all = 10.^([-1.2:0.1:-0.3]);
% rw_all = [0.140:0.001:0.146];
% rw_all = 0.143;

% shear_all = [0:1e-4:shear_Ri0_25]; % Ri=1, shear = 0.97e-3;
% rw_all= 10.^([-2:0.1:-1.1 -1.05 -1:0.0025:-0.7 -0.6:0.1:1]);
% rw_all= 10.^([-2:0.3:-1.1 -1.05 -1:0.0025:-0.9 -0.8:0.1:1]);
% rw_all= 10.^([-1:0.0025:-0.9]);
% rw_all = [0.096:0.001:0.143];
% rw_all=[0.1:0.025:0.3]
% rw_all = [0.142:0.0005:0.145]






    % maxgrow = max(grow)
    % 
    % fig=figure('visible','off');
    % clf;set(gcf,'Color','w');
    % plot(kx_all,grow,'LineWidth',2)
    % xlabel('\it{k_x} (m^{-1})')
    % % title(['Growth rate (1/hour), \it{m_0}=' num2str(m0) ' (m^{-1})'])
    % % plot((rw_all)*N/omega,grow,'LineWidth',2)
    % % xlabel('Wavenumber ratio (k_x/m_z)')
    % xlabel('Natural frequency/Forcing frequency (N k_x/m_z/\omega)')
    % title('Growth rate (1/hour)')
    % grid on;grid minor;
    % ylabel('(1/hour)')
    % set(gca,'fontsize',20)
    % set(gcf, 'InvertHardcopy', 'off')
    % saveas(fig,[[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '_m0' num2str(m0) '.jpeg']);




