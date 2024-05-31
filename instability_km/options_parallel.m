
%%%%% All variables are dimensional variables
clear; close all;
Diffusion = true;
ConvectiveAdjustment = false;
nt_percycle = 72*10; 


%%%%%% exps_flat_diff %%%%%%
expdir = 'experiments_flat_nu2e-4_new/'; %%% Flat bottom with diffusion/viscous dissipation
topo=0;
N = sqrt(1)*1e-3;
Ptide = 43200;
omega = 2*pi/Ptide;
shear_Ri0_25 = 2*N;
shear_Ri1 = N;
shear_all = [0:1e-4:shear_Ri0_25];
m0max = 2*pi/1;
m0min = 2*pi/5000;
k0max = 2*pi/3;
k0min = 2*pi/100000;
m0_all = [0 m0min*[1:1/2:7.5] 0.01:0.01/10:0.6 0.61:0.01/2:1];
kx_all = [0 k0min*[1:1/2:15] 0.001:0.001/10:0.06 0.061:0.001/2:0.1];
lam_z_all = 2*pi./m0_all;
lam_x_all = 2*pi./kx_all;


mkdir(expdir);

% for ns =1:length(shear_all)
parfor ns =11:12
    ns
    % rw_all = rw_mg(ns)
    shear = shear_all(ns)
    
    mkdir([expdir 'shear_' num2str(shear*1e3,3)]);

    rs = shear/omega; %%% shear over omega 
    if(omega==0)
        rs = 0;
    end
   
    for m=1:length(m0_all)
        m
	    m0 = m0_all(m);

    for i=1:length(kx_all)
        kx=kx_all(i);

    % for i=1:length(rw_all)
        % rw = rw_all(i);
        % kx = m0*rw;

        NTtide = 15;
        if(omega==0)
            NTtide = 1/rw/shear/Ptide*10;
        end
        constants;
        loop;
        if(grow(i)>0)
            NTtide = 100;
            constants;
            loop;
        end
        if(grow(i)>0 && grow(i)<0.1)
            NTtide = 400;
            constants;
            loop;
        end
    end

    %%% Save the data
    clear fig a1_t angle_front ct fit_span mz_t pe st tt xx_plot yy_plot buoy dbdt dzetadt ke kew psi re_buoy re_uuu re_www uuu www zeta ke_nond ps_nond
    save([[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '_m0' num2str(m0) '.mat'])

    end

    % plot_timeseires
    % saveas(fig,[expdir '/figs/timeseriesN' num2str(ns) '.jpeg']);
    

end


