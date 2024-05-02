

    prodir = '/Users/ysi/MITgcm_shear_convec/products/';
    expdir = '/Users/ysi/MITgcm_shear_convec/exps_flat'
    % expdir = '/Users/ysi/MITgcm_shear_convec/exps_topo4'


    EXPNAME = {...
        % 'lores_topo0_s0.0015_dz3dx15n-10_lfac1'
        % 'lores_topo0_s0.0015_dz3dx15n-10rand_lfac1'
        % 'lores_topo0_s0.002_dz3dx15n-10rand_lfac1'
        ...
% 'lores_topo4_s0_dz3dx15n-20'
% 'lores_topo4_s0.0001_dz3dx15n-20'
% 'lores_topo4_s0.0002_dz3dx15n-20'
% 'lores_topo4_s0.0003_dz3dx15n-20'
% 'lores_topo4_s0.0004_dz3dx15n-20'
% 'lores_topo4_s0.0005_dz3dx15n-20'
% 'lores_topo4_s0.0006_dz3dx15n-20'
% 'lores_topo4_s0.0007_dz3dx15n-20'
% 'lores_topo4_s0.0008_dz3dx15n-20'
% 'lores_topo4_s0.0009_dz3dx15n-20'
% 'lores_topo4_s0.001_dz3dx15n-20'
% 'lores_topo4_s0.0011_dz3dx15n-20'%%% No grow
% 'lores_topo4_s0.0012_dz3dx15n-20'
% 'lores_topo4_s0.0013_dz3dx15n-20'
% 'lores_topo4_s0.0014_dz3dx15n-20'
% 'lores_topo4_s0.0015_dz3dx15n-20'%%% No grow
% 'lores_topo4_s0.0016_dz3dx15n-20'
% 'lores_topo4_s0.0017_dz3dx15n-20'%%% No grow
% 'lores_topo4_s0.0018_dz3dx15n-20'%%% No grow
% 'lores_topo4_s0.0019_dz3dx15n-20'%%% No grow
% 'lores_topo4_s0.002_dz3dx15n-20'
...
'lores_topo0_s0_dz3dx15n-20'
'lores_topo0_s0.0001_dz3dx15n-20'
'lores_topo0_s0.0002_dz3dx15n-20'
'lores_topo0_s0.0003_dz3dx15n-20'
'lores_topo0_s0.0004_dz3dx15n-20'
'lores_topo0_s0.0005_dz3dx15n-20'
'lores_topo0_s0.0006_dz3dx15n-20'
'lores_topo0_s0.0007_dz3dx15n-20'
'lores_topo0_s0.0008_dz3dx15n-20'
'lores_topo0_s0.0009_dz3dx15n-20'
'lores_topo0_s0.001_dz3dx15n-20' %%% No grow
'lores_topo0_s0.0011_dz3dx15n-20'%%%
'lores_topo0_s0.0012_dz3dx15n-20'%%%
'lores_topo0_s0.0013_dz3dx15n-20'%%%
'lores_topo0_s0.0014_dz3dx15n-20'%%%
'lores_topo0_s0.0015_dz3dx15n-20'%%%
'lores_topo0_s0.0016_dz3dx15n-20'
'lores_topo0_s0.0017_dz3dx15n-20'%%%
'lores_topo0_s0.0018_dz3dx15n-20'%%%
'lores_topo0_s0.0019_dz3dx15n-20'
'lores_topo0_s0.002_dz3dx15n-20'%%%
        };
    nEXP = length(EXPNAME);

