clear;
% expdir = 'eigenvalues/eig_topo4_kv2e-4/';
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

lam_z_all = [1:1:500];
m0_all = 2*pi./lam_z_all;

lam_x_all = [5 10:10:12000];
kx_all = 2*pi./lam_x_all;

Nk = length(kx_all);
Nm = length(m0_all);



m0overk0 = NaN.*zeros(Nk,Nm);
for i=1:Nk
    for j=1:Nm
        m0overk0(i,j) = m0_all(j)/kx_all(i);
    end
end

omega = 2*pi/43200;


for ns=1:Ns
% for ns=11
    shear = shear_all(ns);
    % load([expdir 'shear' num2str(shear*1e3,3) '_output.mat'],'grow','lambda');
    load([expdir 'shear' num2str(shear*1e3,3) '_output.mat'],'grow');

    shear_over_omega = shear/omega;

    exclude_negative_mt = ones(Nk,Nm);
    %%%% Attemp1
    % exclude_negative_mt(m0overk0<shear_over_omega)=NaN;
    %%%% Attempt2
    % exclude_negative_mt(:,1:6)=NaN;
    % exclude_negative_mt(1:2,:)=NaN;
    %%%% Attempt3
    % exclude_negative_mt(:,200:500)=NaN;
    exclude_negative_mt(32:end,:)=NaN;
    %%%% Attempt4
    % Lmax = 250;
    % Lmin = 5;
    % for i=1:Nk
    %     k0=kx_all(i);
    %     for j=1:Nm
    %         m0=m0_all(j);
    %         if((m0-k0*shear/omega)<2*pi/Lmax)
    %             exclude_negative_mt(i,j)=NaN;
    %         end
    %         if((m0-k0*shear/omega)>2*pi/Lmin)
    %             exclude_negative_mt(i,j)=NaN;
    %         end
    %     end
    % end
   
    grow_exclude = exclude_negative_mt.*grow;
    % lambda(:,:,1) = exclude_negative_mt.*lambda(:,:,1);
    % lambda(:,:,2) = exclude_negative_mt.*lambda(:,:,2);

    grow = log(grow)/43200*3600;

    grow_exclude = log(grow_exclude)/43200*3600;

    max_grow_floquet(ns) =max(real(grow),[],'all');
    max_grow_floquet_exclude(ns) =max(real(grow_exclude),[],'all');

    % figure(2);set(gcf,'Color','w')
    % pcolor(kx_all,m0_all,real(grow)');
    % shading flat;colorbar;colormap(redblue);
    % clim([-1 1]*0.3);set(gca,'fontsize',17)
    % xlim([0 0.06]);ylim([0 0.5])
    % xlabel('k_0 (1/m)')
    % ylabel('m_0 (1/m)')
    % % title('The real part of the eigenvalues')
    % title('Growth rate (1/hour)')
    % 
    % figure(3);set(gcf,'Color','w')
    % pcolor(lam_x_all,lam_z_all,real(grow)');
    % shading flat;colorbar;colormap(redblue);
    % clim([-1 1]*0.3);set(gca,'fontsize',17)
    % xlim([0 3000]);ylim([0 250])
    % xlabel('Lx (1/m)')
    % ylabel('Lz (m)')
    % % title('The real part of the eigenvalues')
    % title('Growth rate (1/hour)')

    
    % lambda = real(lambda);
    % grow1 = max(lambda(:,:,1),[],'all');
    % grow2 = max(lambda(:,:,2),[],'all');
    % grow1 = log(grow1)/43200*3600;
    % grow2 = log(grow2)/43200*3600;
    % 
    % grow_floquet1(ns) =max(grow1,[],'all');
    % grow_floquet2(ns) =max(grow2,[],'all');

end



%%
figure(1)
clf;
set(gcf,'Color','w')
plot(shear_all,max_grow_floquet,'LineWidth',2)
hold on;
plot(shear_all,max_grow_floquet_exclude,'LineWidth',2)
% plot(shear_all,grow_floquet1,'--','LineWidth',2)
% % plot(shear_all,grow_floquet2,'LineWidth',2)
set(gca,'fontsize',17)
xlabel('shear (1/s)')
ylabel('(1/hour)')
title('The Floquet exponents as a function of shear')
ylim([0 0.3])
grid on;grid minor;

save('../figures/fig4/Floquet_km_flat.mat','max_grow_floquet','shear_all','max_grow_floquet_exclude')



