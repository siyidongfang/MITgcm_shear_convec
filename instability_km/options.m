
%%%%% All variables are dimensional variables

clear; close all;

Diffusion = false;
ConvectiveAdjustment = false;
% nt_percycle = 72*50; 
nt_percycle = 10;

expdir = 'exps_KHinstability/';
mkdir(expdir);

%%% exps_KHinstability
topo = 0;
omega = 0;
N = sqrt(1)*1e-3;
Ptide = 43200; 
shear_Ri0_25 = 2*N;

%%% exps_Radko19
% topo = 0;
% N = sqrt(10)*1e-3;
% omega = 0.1*1e-3;
% Ptide = 2*pi/omega;
% shear_Ri0_25 = 0.0050; %%% calculated by the script calc_shear_from_Ri

%%% Our simulation
% topo=4;
% N = sqrt(1)*1e-3;
% Ptide = 43200;
% omega = 2*pi/Ptide;
% shear_Ri0_25 = 0.0018;
% % load('rw_mg.mat') %%% The wavenumber ratio m0/k that corresponds to the largest growth rate with out diffusion; or load('rw_mg_test3.mat') with diffusion kappa=1e-5,nu=1e-6.

% shear_all = [0:4e-4:2*shear_Ri0_25]; % Ri=1, shear = 0.97e-3;

shear_all = 0.8*shear_Ri0_25;
 
% rw_all= 10.^([-2:0.1:-1 -0.95:0.01:-0.5 0.6:0.1:1]);
% rw_all= 10.^([-2:0.1:-1.1 -1:0.005:-0.7 -0.6:0.1:1]);
% rw_all = 10.^([-1.2:0.01:-0.5]);
% rw_all = 10.^([-1.2:0.1:-0.3]);
% rw_all = [0.140:0.001:0.146];
% rw_all = 0.143;

rw_all = 10.^([-5:0.25:3]); %%% For K-H instability

m0 =1;

for ns =1:length(shear_all)
    ns
    % rw_all = rw_mg(ns)

    shear = shear_all(ns)
    
    rs = shear/omega; %%% shear over omega 
    if(omega==0)
        rs = 0;
    end
   
    % for i=1:length(rw_all)
    for i = 1
        rw = rw_all(i);
        kx = m0*rw;
        NTtide = 200000;
        constants;
        loop;
        if(grow(i)>0)
            NTtide = 100;
            constants;
            loop;
        end
        % if(grow(i)>0 && grow(i)<1e-2)
        %     NTtide = 30;
        %     constants;
        %     loop;
        % end
        % if(grow(i)>0 && grow(i)<1e-2)
        %     NTtide = 50;
        %     constants;
        %     loop;
        % end
        % if(grow(i)>0 && grow(i)<1e-2)
        %     NTtide = 50;
        %     constants;
        %     loop;
        % end

    end

    % plot_timeseires
    
    %%% Save the data
    clear buoy dbdt dzetadt ke kew psi re_buoy re_uuu re_www uuu www zeta
    save([expdir '/growth_shear' num2str(shear*1e3,3) '.mat'])

    maxgrow = max(grow)
    
    figure(1)
    clf;set(gcf,'Color','w');
    semilogx((rw_all),grow,'LineWidth',2)
    grid on;grid minor;
    title('Growth rate (1/hour)')
    ylabel('(1/hour)')
    xlabel('Wavenumber ratio (k_x/m_z)')
    set(gca,'fontsize',20)
    % set(gcf, 'InvertHardcopy', 'off')
    print('-djpeg','-r150',[expdir '/growth_shear' num2str(shear*1e3,3) '.jpeg']);

end




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





