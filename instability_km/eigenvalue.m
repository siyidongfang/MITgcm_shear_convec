% Define the symbolic variable t
% syms t sigma a1(t) a2(t) a0 nu kappa

clear;close all;
expdir = 'eigenvalues/';

nu=2e-4;
kappa=2e-4;
Ptide = 43200;
% pi = sym(pi);
omega = 2*pi/Ptide;
topo=4;
cs = cosd(topo);
ss = sind(topo);
N = 1e-3;

topo=4;
shear_Ri0_25 = 0.0018; % 0.0017525;
shear_Ri1 = 9.7e-04;

% topo=0;
% shear_Ri0_25 = 2*N;
% shear_Ri1 = N;

shear_all = [0:1e-4:shear_Ri0_25]; 

lz_all = [1:1:500];
m0_all = 2*pi./lz_all;

lx_all = [5 10:10:12000];
k0_all = 2*pi./lx_all;

NTtide = 1;
nt_percycle = 36;
dt = Ptide/nt_percycle;
Nt = round(NTtide*Ptide/dt);
tt = dt:dt:Nt*dt;

Nk = length(k0_all);
Nm = length(m0_all);
Ns = length(shear_all);

real_eigv = NaN*ones(Nk,Nm,Nt);

for ns = 1:Ns
    shear = shear_all(ns);
    rs = shear/omega;
    parfor i=1:Nk
        % i
        k0=k0_all(i);
    
        for j=1:Nm
            m0 = m0_all(j);
        
            for o=1:Nt
                t0 =tt(o);
                st = sin(omega*t0);
                
                mz = m0-rs*st*k0;
                a0 = (1i*m0*ss-1i*k0*cs)*N^2; 
                a1 = -(k0^2+mz^2); 
                a2 = 1i*k0*cs - 1i*mz*ss;
                
                % Define the matrix M(t)
                M = [nu*a1, a2; a0/a1, kappa*a1];
                % Compute the eigenvalues and eigenvectors
                [V, D] = eig(M);
                % The eigenvalues are the diagonal elements of D
                sigma = diag(D);
    
                real_eigv(i,j,o) = real(max(sigma));
            
            end
        
        end
    end
    

  
    sigma_all = real_eigv;
    sigma_all_mean = mean(sigma_all,3)*3600;
    sigma_all_max = max(sigma_all,[],3)*3600;

    figure(1)
    % pcolor(k0_all,m0_all,sigma_all_max');shading flat;
    pcolor(lx_all/1000,lz_all,sigma_all_max');shading flat;
    colorbar;colormap(WhiteBlueGreenYellowRed(0))
    % clim([0 1]*3e-5)
    clim([0 1]*0.4)
    % xlim([0 1])
    % ylim([0 1])
    set(gca,'Fontsize',17);
    xlabel('Across-slope wavelength (km)')
    ylabel('Slope-normal wavelength (m)')
    title('Maximum Real\{\sigma\} in one tidal cycle (1/hour)')
    print('-dpng','-r200',[expdir 'topo4_shear' num2str(shear*1e3,3) '.png']);

    
    max_growth = max(sigma_all_max,[],'all');
    
    % figure(1)
    % aaa = squeeze(real_eigv(19,300,:))
    % plot(tt,aaa)

    save([expdir 'topo4_shear' num2str(shear*1e3,3) '_output.mat']);

end



