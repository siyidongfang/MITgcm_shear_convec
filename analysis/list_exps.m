

    prodir = '/Users/ysi/MITgcm_shear_convec/products/';
    % expdir = '/Users/ysi/MITgcm_shear_convec/exps_flat/';
    % expdir = '/Users/ysi/MITgcm_shear_convec/exps_flat_hires/';
    expdir = '/Users/ysi/MITgcm_shear_convec/exps_topo4_hires/';


    EXPNAME = {...
        % 'H1500_lfac2_smoohalf120m_lores_topo0_s0.0006_dz3dx20n-20'
        % 'H1500_lores_topo0_s0.0006_dz3dx20n-20'
        % 'H1500_smoohalf120m_lores_topo0_s0.0006_dz3dx20n-20'
        % ...
        % 'lores_topo0_s0.0005_dz3dx15n-20'
        % 'lores_topo0_s0.0006_dz3dx15n-20rand_lfac1_smooth'
        ...
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
        'hires_topo4_s0_dz1dx6n-20'
        'hires_topo4_s0.0001_dz1dx6n-20'
        'hires_topo4_s0.0002_dz1dx6n-20'
        'hires_topo4_s0.0003_dz1dx6n-20'
        'hires_topo4_s0.0004_dz1dx6n-20'
        'hires_topo4_s0.0005_dz1dx6n-20'
        'hires_topo4_s0.0006_dz1dx6n-20'
        'hires_topo4_s0.0007_dz1dx6n-20'
        'hires_topo4_s0.0008_dz1dx6n-20'
        'hires_topo4_s0.0009_dz1dx6n-20'
        'hires_topo4_s0.001_dz1dx6n-20'
        'hires_topo4_s0.0011_dz1dx6n-20'
        'hires_topo4_s0.0012_dz1dx6n-20'
        'hires_topo4_s0.0013_dz1dx6n-20'
        'hires_topo4_s0.0014_dz1dx6n-20'
        'hires_topo4_s0.0015_dz1dx6n-20'
        'hires_topo4_s0.0016_dz1dx6n-20'
        'hires_topo4_s0.0017_dz1dx6n-20'
        'hires_topo4_s0.0018_dz1dx6n-20'
        'hires_topo4_s0.0019_dz1dx6n-20'
        'hires_topo4_s0.002_dz1dx6n-20'
        };
    nEXP = length(EXPNAME);

