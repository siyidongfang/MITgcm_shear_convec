clear;
% expdir = 'eigenvalues/eig_topo4_kv2e-4/';
% expdir = 'eigenvalues/eig_flat_kv2e-4-lores/';
expdir = 'eigenvalues/eig_flat_kv2e-4/';

N = sqrt(1)*1e-3;
% topo=4;
% shear_Ri0_25 = 0.0018; % 0.0017525;
% shear_Ri1 = 9.7e-04;

topo=0;
shear_Ri0_25 = 2*N;
shear_Ri1 = N;

shear_all = [0:1e-4:shear_Ri0_25]; 
Ns = length(shear_all);

% lam_z_all = [1:1:500];
lam_z_all = [0.01:0.01:0.09 0.1:0.1:10 10.25:0.25:50 50.5:0.5:500];
m0_all = 2*pi./lam_z_all;

% lam_x_all = [5 10:10:12000];
lam_x_all = [0.05 0.1 0.25:0.25:10 10:1:800 805:5:12000];
kx_all = 2*pi./lam_x_all;


Nk = length(kx_all);
Nm = length(m0_all);

xxx = [1:1:5000];
yyy924 = xxx/9.24;
yyy9 = xxx/9;
yyy696 = xxx/6.96;




% for ns=1:Ns
for ns=11
    shear = shear_all(ns);
    load([expdir 'shear' num2str(shear*1e3,3) '_output.mat'],'grow');
    grow = real(log(grow)/43200*3600);
    max_grow_floquet(ns) =max(real(grow),[],'all');

    grow(grow<0)=NaN;
    % figure(2);set(gcf,'Color','w')
    % pcolor(kx_all,m0_all,grow');
    % shading flat;colorbar;colormap(redblue);
    % clim([-1 1]*0.3);set(gca,'fontsize',17)
    % xlim([0 0.06]);ylim([0 0.5])
    % xlabel('k_0 (1/m)')
    % ylabel('m_0 (1/m)')
    % % title('The real part of the eigenvalues')
    % title('Growth rate (1/hour)')

    figure(3);
    clf
    set(gcf,'Color','w')
    pcolor(lam_x_all,lam_z_all,grow');
    hold on;
    plot(xxx,yyy924,'--','LineWidth',2,'Color','k')
    plot(xxx,yyy9,'--','LineWidth',2,'Color','k')
    plot(xxx,yyy696,'--','LineWidth',2,'Color','k')
    shading flat;colorbar;
    % colormap(redblue);clim([-1 1]*0.1);
    colormap(jet);clim([0 1]*0.06);
    set(gca,'fontsize',17)
    xlim([0 6000]);
    % ylim([0 250])
    xlabel('Lx (1/m)')
    ylabel('Lz (m)')
    % title('The real part of the eigenvalues')
    title('Growth rate (1/hour)')
    grid on;grid minor;

end


