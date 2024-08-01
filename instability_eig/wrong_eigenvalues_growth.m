
clear;close all;
expdir = 'eigenvalues/eigenvalues_wrong/';

N = 1e-3;

topo=4;
shear_Ri0_25 = 0.0018; % 0.0017525;
shear_Ri1 = 9.7e-04;

% topo=0;
% shear_Ri0_25 = 2*N;
% shear_Ri1 = N;

shear_all = [0:1e-4:shear_Ri0_25]; 
Ns = length(shear_all);

lz_all = [1:1:500];
m0_all = 2*pi./lz_all;

lx_all = [5 10:10:12000];
k0_all = 2*pi./lx_all;
Nk = length(k0_all);
Nm = length(m0_all);

zeta = 1e-20*ones(Nk,Nm);

for ns = 1:Ns
    ns
    shear = shear_all(ns);
    load([expdir 'topo4_shear' num2str(shear*1e3,3) '_output.mat'],'max_growth','sigma_all','dt','Nt');
    grow_eigv(ns) = max_growth;

    % parfor i=1:Nk
    %     % i
    %     k0=k0_all(i);
    % 
    %     for j=1:Nm
    %         m0 = m0_all(j);
    %         for o=2:Nt
    %             zeta(i,j) = zeta(i,j) + sigma_all(i,j,o-1)*zeta(i,j)*dt;
    %         end
    %     end
    % end
    % save([expdir 'topo4_shear' num2str(shear*1e3,3) '_zeta.mat'],'zeta');
end



%%
% for ns = 1:Ns
%     shear = shear_all(ns);
%     load([expdir 'topo4_shear' num2str(shear*1e3,3) '_zeta.mat'],'zeta');
%     max_zeta(ns) = max(zeta,[],'all');
% end
% 
% growth_zeta(ns) =log10(max_zeta(ns));

figure(2)
set(gcf,'Color','w')
% plot(shear_all,growth_zeta)
plot(shear_all,grow_eigv,'LineWidth',2)
set(gca,'fontsize',17)
xlabel('shear (1/s)')
ylabel('(1/hour)')
title('Maximum Real\{\sigma (t,k_0,m_0)\} as a function of shear')
ylim([0 0.3])
grid on;grid minor;



