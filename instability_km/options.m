
%%%%% All variables are dimensional variables

clear; close all;
Diffusion = true;
ConvectiveAdjustment = false;
nt_percycle = 720; 

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


%%%%%% exps_Radko %%%%%%
expdir = 'exps_Radko/';
topo = 0;
N = sqrt(10)*1e-3;
omega = 0.1*1e-3;
Ptide = 2*pi/omega;
shear_Ri0_25 = 0.0050; %%% calculated by the script calc_shear_from_Ri
shear_Ri1 = 0.0031625;
shear_all = shear_Ri1;
m0_all = [0:0.01:6];
kx_all = [-0.1:0.00025:0.1];

% m0_all = [0:3:6];
% kx_all = [-0.1:0.05:0.1];

% %%%%%% exps_test %%%%%%
% expdir = 'exps_test/';
% topo=4;
% N = sqrt(1)*1e-3;
% Ptide = 43200;
% omega = 2*pi/Ptide;
% shear_Ri0_25 = 0.0018;
% shear_Ri1 = 0.97e-3;
% % load('rw_mg.mat') %%% The wavenumber ratio m0/k that corresponds to the largest growth rate with out diffusion; or load('rw_mg_test3.mat') with diffusion kappa=1e-5,nu=1e-6.
% shear_all = [0:0.1e-4:shear_Ri0_25]; % Ri=1, shear = 0.97e-3;
% rw_all= 10.^([-2:0.1:-1.1 -1.05 -1:0.0025:-0.7 -0.6:0.1:1]);
% m0 =1;

mkdir(expdir);

for ns =1:length(shear_all)
    ns
    % rw_all = rw_mg(ns)
    shear = shear_all(ns)
    
    rs = shear/omega; %%% shear over omega 
    if(omega==0)
        rs = 0;
    end
   
    for m=1:length(m0_all)
            m0 = m0_all(m);


        for i=1:length(kx_all)
            kx=kx_all(i);

        % for i=1:length(rw_all)
        %     rw = rw_all(i);
            % kx = m0*rw;
    
            NTtide = 10;
            if(omega==0)
                NTtide = 1/rw/shear/Ptide*10;
            end
            constants;
            loop;
            if(grow(i)>0)
                NTtide = 30;
                constants;
                loop;
            end
            if(grow(i)>0 && grow(i)<5e-2)
                NTtide = 100;
                constants;
                loop;
            end
            if(grow(i)>0 && grow(i)<5e-3)
                NTtide = 600;
                constants;
                loop;
            end
            if(grow(i)>0 && grow(i)<1e-4)
                NTtide = 1000;
                constants;
                loop;
            end
        end


    %%% Save the data
    clear a1_t angle_front ct fit_span mz_t pe st tt xx_plot yy_plot buoy dbdt dzetadt ke kew psi re_buoy re_uuu re_www uuu www zeta ke_nond ps_nond
    save([expdir '/growth_shear' num2str(shear*1e3,3) '_m0' num2str(m0) '.mat'])

    maxgrow = max(grow)
    
    fig=figure(1);
    set(fig,'visible','off');
    clf;set(gcf,'Color','w');
    plot(kx_all,grow,'LineWidth',2)
    xlabel('\it{k_x} (m^{-1})')
    title(['Growth rate (1/hour), \it{m_0}=' num2str(m0) ' (m^{-1})'])
    % semilogx((rw_all),grow,'LineWidth',2)
    % xlabel('Wavenumber ratio (k_x/m_z)')
    % title('Growth rate (1/hour)')
    grid on;grid minor;
    ylabel('(1/hour)')
    set(gca,'fontsize',20)
    set(gcf, 'InvertHardcopy', 'off')
    saveas(fig,[expdir '/growth_shear' num2str(shear*1e3,3) '_m0' num2str(m0) '.jpeg']);


    end

    % plot_timeseires


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
% save(['output_topo0/growth_topo' num2str(topo) '_mz' num2str(mz) 'kx' num2str(kx) '_NOdiff.mat'])



% rw_all= 10.^([-2:0.1:-1 -0.95:0.01:-0.5 0.6:0.1:1]);
% rw_all= 10.^([-2:0.1:-1.1 -1:0.005:-0.7 -0.6:0.1:1]);
% rw_all = 10.^([-1.2:0.01:-0.5]);
% rw_all = 10.^([-1.2:0.1:-0.3]);
% rw_all = [0.140:0.001:0.146];
% rw_all = 0.143;



