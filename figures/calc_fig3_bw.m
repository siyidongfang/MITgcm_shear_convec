
clear;close all;

%--- Find out the shear corresponding to this Ri
Ri_set = 0.6;
load('fig3/Ri_topo4.mat','Ri_min','shear_calc_Ri')
[a b] = min(abs(Ri_set - Ri_min));
Ri = Ri_min(b);
shear_set = shear_calc_Ri(b);

%--- load one simulation with Ri=0.4 (Inverse Ri = 2.5)
load('fig3/topo4_kappa0.mat','shear_all','crop_limit')
[a num_shear] = min(abs(shear_set - shear_all));
shear = shear_all(num_shear);
expdir = '../instability_km/exps/topo4_kappa0/';

%--- find out the most unstable mode
h_shear = 2000;
m0_rw = 2*pi/h_shear;
m0max = 2*pi/1;
m0min = m0_rw;
lam_x_all = [10:10:100 125:25:1000 1100:100:5000 5200:200:10000 10000:500:40000 40000:1000:100000];
lam_x_all = flip(lam_x_all);
kx_all = 2*pi./lam_x_all;
rw_all = kx_all/m0_rw;
lam_x_real = 2*pi./(m0_rw.*rw_all);

Nrw = length(rw_all);
for i=1:Nrw
    fname = [[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '_rw' num2str(i) '.mat'];
    load(fname,'grow');
    grow_r(i)=grow(i);
end


rw_idx=find(lam_x_real<=crop_limit);

grow_r_crop = grow_r(rw_idx);
lam_x_real_crop = lam_x_real(rw_idx);
rw_all_crop = rw_all(rw_idx);

figure(1)
plot(lam_x_real,grow_r)
hold on;
plot(lam_x_real(rw_idx),grow_r(rw_idx))

[a b] = max(grow_r_crop);
num_rw = rw_idx(b)
max_grow = grow_r_crop(b)
rw_m = rw_all_crop(b)
kx_m = rw_m*m0_rw
lam_x_m = 2*pi./kx_m




