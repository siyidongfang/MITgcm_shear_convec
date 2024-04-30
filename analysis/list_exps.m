

    prodir = '/Users/ysi/MITgcm_shear_convec/products/';
    expdir = '/Users/ysi/MITgcm_shear_convec/exps_test'

    EXPNAME = {...
        'noBtide_nosmooth_topo4_s0_dz3dx15n1e-6'
        'noBtide_nosmooth_topo4_s0.0001_dz3dx15n1e-6'
        'noBtide_nosmooth_topo4_s0.0015_dz3dx15n1e-6'
        'fixTide_smooth_topo4_s0_dz3dx15n1e-15_lfac1'
        'fixTide_smooth_topo4_s0.0001_dz3dx15n1e-15'
        'fixTide_smooth_topo4_s0.0005_dz3dx15n1e-15'
        'fixTide_smooth_topo4_s0_dz3dx15n1e-15'
        'nosmooth_topo4_s0_dz3dx15n1e-6'
        'nosmooth_topo4_s0.0001_dz3dx15n1e-6'
        'nosmooth_topo4_s0.0015_dz3dx15n1e-6'
        'LinearN_noNoise_lfac1_topo4_s0.0015_dz3dx15n1e-10'
        'LinearN_noNoise_lfac1_topo4_s0_dz3dx15n1e-10'
        'linearN_smooth_topo4_s0_dz3dx15n1e-15'
        };
    nEXP = length(EXPNAME);

