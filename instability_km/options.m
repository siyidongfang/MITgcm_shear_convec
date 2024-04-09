
%%%%% All variables are dimensional variables

clear; close all;

Diffusion = true;
nt_percycle = 72*10; 

% topo = 0;
% N = sqrt(10)*1e-3;
% omega = 0.1*1e-3;
% Ptide = 2*pi/omega;

expdir = 'output/topo4_Nsq1e-6_test3';
mkdir(expdir);

topo=4;
N = sqrt(1)*1e-3;
Ptide = 43200;
omega = 2*pi/Ptide;

% calc_shear_from_Ri
% shear_Ri0_25 = 0.0050;
shear_Ri0_25 = 0.0018;
shear_all = [0:1e-4:shear_Ri0_25]; % Ri=1, shear = 0.97e-3;

% rw_all= 10.^([-2:0.1:-1 -0.95:0.01:-0.5 0.6:0.1:1]);
% rw_all= 10.^([-2:0.1:-1.1 -1:0.005:-0.7 -0.6:0.1:1]);
% rw_all = 10.^([-1.2:0.01:-0.5]);
% rw_all = 10.^([-1.2:0.1:-0.3]);
% rw_all = [0.140:0.001:0.146];
% rw_all = 0.143;
m0 =1;
load('rw_mg_test3.mat')


for ns =1:length(shear_all)
    ns
    rw_all = rw_mg(ns)

    shear = shear_all(ns)
    rs = shear/omega; %%% shear over omega 
   
    for i=1:length(rw_all)
        i
        rw = rw_all(i);
        kx = m0*rw;
        NTtide = 200;
        constants;
        loop;
        % if(grow(i)>0 && grow(i)<1e-3)
        % if(grow(i)>0)
        %     NTtide = 150;
        %     constants;
        %     loop;
        % end
        % if(grow(i)>0)
        %     NTtide = 400;
        %     constants;
        %     loop;
        % end
        % if(grow(i)>0)
        %     NTtide = 600;
        %     constants;
        %     loop;
        % end
    end
    
    %%% Save the data
    % clear buoy dbdt dzetadt ke kew psi re_buoy re_uuu re_www uuu www zeta
    save([expdir '/growth_shear' num2str(shear*1e3,3) '_test3_analysis.mat'])
    
    % figure(1)
    % clf;set(gcf,'Color','w');
    % semilogx((rw_all),grow,'LineWidth',2)
    % grid on;grid minor;
    % title('Growth rate (1/hour)')
    % ylabel('(1/hour)')
    % xlabel('Wavenumber ratio log10(k_x/m_z)')
    % set(gca,'fontsize',20)
    % % set(gcf, 'InvertHardcopy', 'off')
    % print('-djpeg','-r150',[expdir '/growth_shear' num2str(shear*1e3,3) '.jpeg']);

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





