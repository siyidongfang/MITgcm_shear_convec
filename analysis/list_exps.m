

    expdir = '/Users/ysi/MITgcm_shear_convec/exps_test/';
    prodir = '/Users/ysi/MITgcm_shear_convec/products/';
    EXPNAME = {...
        'shear6e-4Hs300_Nr241Nx100_spin'
        'shear7e-4Hs300_Nr241Nx100_spin'
        'shear8e-4Hs300_Nr241Nx100_spin'
        'shear9e-4Hs300_Nr241Nx100_spin'
        'shear7e-4Hs300_Nr241Nx100_spin_noise0.01'
        'shear1e-3Hs300_Nr241Nx100_spin-lores'
        'test_bottomGrid'
        'test_halfw_bottom'
        'shear1e-3Hs300_Nr241Nx100_spin_iniN2_1e-10'
        'shear1e-3Hs300_Nr241Nx100_spin_ini0'
        'shear1e-3Hs300_Nr241Nx200_spin'
        'test_bottomGrid-hires'
        'test_addTide2presure_noSurfaceRelax'
        'test_turnOffRestoring'
        'test_CorioInHydroPressure'
        'test_gv_fwsin'
        'test_N2sin2_surfacesponge'
        'test_surfacesponge_v1.3'
        'noNoise_noCadj'
        'noNoise_cadj'
        'noNoise_cadj_diff'
        'noNoise_kpp'
        'noise1e-3_noCadj'
        'noise1e-3_cadj'
        'noise1e-3_cadj_diff'
        'noise1e-3_kpp'
        'ushear3e-4Hs200_Nr376Nx320_spin'
        };
    nEXP = length(EXPNAME);
