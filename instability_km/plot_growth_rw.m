
clear;
load('output_Ri1_Nsq1e-6_topo4/growth_topo4_NOdiff.mat')
figure(1)
semilogx(rw_all,growth)

[maxgrow rwidx] = max(growth)

rw_maxgrow = rw_all(rwidx);
kx = mz*rw_maxgrow;

load('/Users/ysi/MITgcm_shear_convec/instability_km/output_Ri1_Nsq1e-6_topo4/mz0.01kx0.0015849shear0.00097.mat')
plot_timeseires


