

    prodir = '/Users/ysi/MITgcm_shear_convec/products/';
    expdir = '/Users/ysi/MITgcm_shear_convec/exps_flat_lores/';
    % expdir = '/Users/ysi/MITgcm_shear_convec/exps_flat/';
    % expdir = '/Users/ysi/MITgcm_shear_convec/exps_topo4/';


    EXPNAME = { ... 
'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_smag_h1e-4v2e-4'
'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_smag_h1e-4v5e-4'
'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_smag_h1e-4v1e-3'
'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_smag_h5e-4v2e-4'
...
'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_nosmag_h1e-4v2e-4'
'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_nosmag_hv2e-4'
'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_nosmag_h2e-4v1e-4'
'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_nosmag_h4e-4v2e-4'
'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_nosmag_h6e-4v3e-4'
'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_nosmag_h1e-3v1e-4'

% 'topo0_H500_smo100m_s0_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0001_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0002_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0003_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0004_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0005_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0006_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0007_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0008_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0009_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.001_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0011_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0012_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0013_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0014_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0015_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0016_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0017_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0018_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.0019_dz1dx3ln200n-20'
% 'topo0_H500_smo100m_s0.002_dz1dx3ln200n-20'
...
% 'topo4_H500_smo100m_s0_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0001_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0002_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0003_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0004_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0005_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0006_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0007_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0008_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0009_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.001_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0011_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0012_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0013_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0014_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0015_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0016_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0017_dz1dx3ln200n-20'
% 'topo4_H500_smo100m_s0.0018_dz1dx3ln200n-20'
        };
    nEXP = length(EXPNAME);



    % expdir = '/Users/ysi/MITgcm_shear_convec/exps_flat_hires/';
    % expdir = '/Users/ysi/MITgcm_shear_convec/exps_topo4_hires/';
    % expdir = '/Volumes/MIT/MITgcm_shear_convec/exps_backup_2024/backup_wrongshear/';


        % 'topo0_H900_smo90m_s0.0006_dz3dx20ln400n-7'
        % 'H900_smooth120m_topo0_s0.0006_dz3dx20n-5'
        % 'H900_smooth120m_topo0_s0.0006_dz3dx20n-6'
        % 'H900_smooth120m_topo0_s0.0006_dz3dx20n-7'
        % 'H950_whitenoise_fixSNX_smooth100m_topo0_s0.0006_dz1dx3n-7'
        % 'H1500_whitenoise_fixSNX_smooth120m_topo0_s0.0006_dz3dx20n-4'
        % 'H1500_whitenoise_fixSNX_smooth120m_topo0_s0.0006_dz3dx20n-9'
        % 'H1500_fixSNX_smooth120m_topo0_s0.0006_dz3dx20n-9'
        % 'H1500_lores_topo0_s0.0006_dz3dx20n-9'
        % 'H1500_lores_topo0_s0.0006_dz3dx20n-9_lessDiag' %%% smooth120m
        % 'H1500_lfac2_smoohalf120m_lores_topo0_s0.0006_dz3dx20n-20'
        % 'H1500_smoohalf120m_lores_topo0_s0.0006_dz3dx20n-20'
        % 'H1500_NOsmooth_topo0_s0.0006_dz3dx20n-12'
        % 'H1500_lores_topo0_s0.0006_dz3dx20n-20' %%% smooth45m
        ...
        % 'lores_topo0_s0.0005_dz3dx15n-20'
        % 'lores_topo0_s0.0006_dz3dx15n-20rand_lfac1_smooth'
        % ...
        % 'hires_topo0_s0_dz1dx6n-20'
        % 'hires_topo0_s0.0001_dz1dx6n-20'
        % 'hires_topo0_s0.0002_dz1dx6n-20'
        % 'hires_topo0_s0.0003_dz1dx6n-20'
        % 'hires_topo0_s0.0004_dz1dx6n-20'
        % 'hires_topo0_s0.0005_dz1dx6n-20'
        % 'hires_topo0_s0.0006_dz1dx6n-20'
        % 'hires_topo0_s0.0007_dz1dx6n-20'
        % 'hires_topo0_s0.0008_dz1dx6n-20'
        % 'hires_topo0_s0.0009_dz1dx6n-20'
        % 'hires_topo0_s0.001_dz1dx6n-20'
        % 'hires_topo0_s0.0011_dz1dx6n-20'
        % 'hires_topo0_s0.0012_dz1dx6n-20'
        % 'hires_topo0_s0.0013_dz1dx6n-20'
        % 'hires_topo0_s0.0014_dz1dx6n-20'
        % 'hires_topo0_s0.0015_dz1dx6n-20'
        % 'hires_topo0_s0.0016_dz1dx6n-20'
        % 'hires_topo0_s0.0017_dz1dx6n-20'
        % 'hires_topo0_s0.0018_dz1dx6n-20'
        % 'hires_topo0_s0.0019_dz1dx6n-20'
        % 'hires_topo0_s0.002_dz1dx6n-20'
        ...
        ...
        % 'hires_topo4_s0_dz1dx6n-20'
        % 'hires_topo4_s0.0001_dz1dx6n-20'
        % 'hires_topo4_s0.0002_dz1dx6n-20'
        % 'hires_topo4_s0.0003_dz1dx6n-20'
        % 'hires_topo4_s0.0004_dz1dx6n-20'
        % 'hires_topo4_s0.0005_dz1dx6n-20'
        % 'hires_topo4_s0.0006_dz1dx6n-20'
        % 'hires_topo4_s0.0007_dz1dx6n-20'
        % 'hires_topo4_s0.0008_dz1dx6n-20'
        % 'hires_topo4_s0.0009_dz1dx6n-20'
        % 'hires_topo4_s0.001_dz1dx6n-20'
        % 'hires_topo4_s0.0011_dz1dx6n-20'
        % 'hires_topo4_s0.0012_dz1dx6n-20'
        % 'hires_topo4_s0.0013_dz1dx6n-20'
        % 'hires_topo4_s0.0014_dz1dx6n-20'
        % 'hires_topo4_s0.0015_dz1dx6n-20'
        % 'hires_topo4_s0.0016_dz1dx6n-20'
        % 'hires_topo4_s0.0017_dz1dx6n-20'
        % 'hires_topo4_s0.0018_dz1dx6n-20'
        % 'hires_topo4_s0.0019_dz1dx6n-20'
        % 'hires_topo4_s0.002_dz1dx6n-20'