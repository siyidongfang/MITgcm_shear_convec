
clear;
NTtide = 1;

expdir = 'K1/';

Diffusion = false;
nt_percycle = 72*30; 

Ptide_all = [1 2]*43200;
topo_all = [0:2:20];
N_all = [0 0.1 0.5 1:2:11]*1e-3;

Ptide = 43200*2;
omega = 2*pi/Ptide;

Nn = length(N_all);
Ntopo = length(topo_all);

h_shear = 250;
m0_rw = 2*pi/h_shear;

Ns = 251; % length of shear_all

kappa_const = 1.4e-7;
nu_const = 1e-6;

if(Diffusion)
    kappa = kappa_const;
    nu = nu_const;

    lam_z_all = 0:1:500;
    lam_x_all = [5 10:10:12000];
else 
    kappa = 0;
    nu = 0;

    lam_z_all = h_shear;
    lam_x_all = [-8000:10:8000];
end

m0_all = 2*pi./lam_z_all;
kx_all = 2*pi./lam_x_all;
Nk = length(kx_all);
Nm = length(m0_all);




Ns = 31;
grow_floquet = NaN.*zeros(Ntopo,Nn,Ns);
max_kidx = NaN.*zeros(Ntopo,Nn,Ns);
max_midx = NaN.*zeros(Ntopo,Nn,Ns);

for ntopo = 1:Ntopo
    topo = topo_all(ntopo)
    cs = cosd(topo);
    ss = sind(topo);
    for nn = 1:Nn
        N = N_all(nn);
        clear shear_all;
        shear_all = [0:3*N/30:3*N]; 
        Ns = length(shear_all);
        for ns =1:Ns
            shear = shear_all(ns);
            fname = [expdir 'ptide' num2str(Ptide) '_topo' num2str(topo) '_N' num2str(N*1e3,3) '_shear' num2str(shear*1e3,3) '.mat'];
            load(fname,'grow'); 
            % Find the most unstable mode and the maximum growth rate
            [grow_floquet(ntopo,nn,ns),I] =max(grow,[],'all');
            [max_kidx(ntopo,nn,ns), max_midx(ntopo,nn,ns)] = ind2sub(size(grow), I);
        end
    end
end

save([expdir 'grow_K1.mat']);

