

    % prodir = '../../MITgcm_shear_convec/products/';
    expdir ='../../MITgcm_shear_convec/exps_flat_hires/'
    % expdir ='../../MITgcm_shear_convec/exps_topo4_hires/'
    % expdir ='../../MITgcm_shear_convec/exps_flat/'
    % expdir ='../../MITgcm_shear_convec/exps_topo4/'

    EXPNAME = {
        'topo0_H500Lx3k_s0.0dz1dx3n-20sm100_kv8e-5' %
        'topo0_H500Lx3k_s0.1dz1dx3n-20sm100_kv8e-5' %
        'topo0_H500Lx3k_s0.2dz1dx3n-20sm100_kv8e-5' %80
        'topo0_H500Lx3k_s0.3dz1dx3n-20sm100_kv8e-5' %
        'topo0_H500Lx3k_s0.4dz1dx3n-20sm100_kv8e-5' %63
        'topo0_H500Lx3k_s0.5dz1dx3n-20sm100_kv8e-5' %
        'topo0_H500Lx3k_s0.6dz1dx3n-20sm100_kv8e-5' %42
        'topo0_H500Lx3k_s0.7dz1dx3n-20sm100_kv8e-5' %
        'topo0_H500Lx3k_s0.8dz1dx3n-20sm100_kv8e-5' %32
        'topo0_H500Lx3k_s0.9dz1dx3n-20sm100_kv8e-5' %
        'topo0_H500Lx3k_s1.0dz1dx3n-20sm100_kv8e-5' %26
        'topo0_H500Lx3k_s1.1dz1dx3n-20sm100_kv8e-5' %
        'topo0_H500Lx3k_s1.2dz1dx3n-20sm100_kv8e-5' %21
        'topo0_H500Lx3k_s1.3dz1dx3n-20sm100_kv8e-5' %
        'topo0_H500Lx3k_s1.4dz1dx3n-20sm100_kv8e-5' %18
        'topo0_H500Lx3k_s1.5dz1dx3n-20sm100_kv8e-5' %
        'topo0_H500Lx3k_s1.6dz1dx3n-20sm100_kv8e-5' %16
        'topo0_H500Lx3k_s1.7dz1dx3n-20sm100_kv8e-5' %
        'topo0_H500Lx3k_s1.8dz1dx3n-20sm100_kv8e-5' %14
        'topo0_H500Lx3k_s1.9dz1dx3n-20sm100_kv8e-5' %
        'topo0_H500Lx3k_s2.0dz1dx3n-20sm100_kv8e-5' %13
        'topo0_H500Lx3k_s2.3dz1dx3n-20sm100_kv8e-5' %
        ...
        };
    nEXP = length(EXPNAME);




         % 'topo0_H500Lx3k_s1.0dz1dx3n-20sm100_kv2e-4'        %49
        % 'topo0_H500Lx2.775k_s1.0dz1dx1.85n-20sm100_kv2e-4' %15
        % 'topo0_H500Lx3.08k_s1.0dz1dx3.08n-20sm100_kv2e-4'  %25
        % 'topo0_H500Lx3.48k_s1.0dz1dx3.48n-20sm100_kv2e-4'  %30
        % 'topo0_H500Lx3.08k_s1.0dz1dx3.08n-20sm100_kv1e-4' %52
        % 'topo0_H500Lx3.08k_s1.0dz1dx3.08n-20sm100_kv7.5e-5'  %50
        % 'topo0_H500Lx3.08k_s1.0dz1dx3.08n-20sm100_kv5e-5' %25



    % expdir ='/Volumes/MIT/MITgcm_shear_convec/exps_topo4/'
    % expdir ='/Volumes/MIT/MITgcm_shear_convec/exps_flat_test/'
    % expdir ='/Volumes/MIT/MITgcm_shear_convec/experiments/'
    
    
    % expdir = '../../MITgcm_shear_convec/exps_hires/';
    % expdir = '../../MITgcm_shear_convec/exps_SteadyShear/';
    % expdir = '../../MITgcm_shear_convec/experiments/';
    % expdir = '../../MITgcm_shear_convec/exps_flat_lores/';
    % expdir = '../../MITgcm_shear_convec/exps_flat/';
    % expdir = '../../MITgcm_shear_convec/exps_topo4/';




        ...
        % 'topo4_H500Lx3k_s0.0dz1dx3n-20sm100_kv2e-4' %40
        % 'topo4_H500Lx3k_s0.2dz1dx3n-20sm100_kv2e-4' %40
        % 'topo4_H500Lx3k_s0.4dz1dx3n-20sm100_kv2e-4' %40
        % 'topo4_H500Lx3k_s0.6dz1dx3n-20sm100_kv2e-4' %40
        % 'topo4_H500Lx3k_s0.8dz1dx3n-20sm100_kv2e-4' %32
        % 'topo4_H500Lx3k_s1.0dz1dx3n-20sm100_kv2e-4' %25
        % 'topo4_H500Lx3k_s1.2dz1dx3n-20sm100_kv2e-4' %21
        % 'topo4_H500Lx3k_s1.4dz1dx3n-20sm100_kv2e-4' %18
        % 'topo4_H500Lx3k_s1.6dz1dx3n-20sm100_kv2e-4' %16
        % 'topo4_H500Lx3k_s1.8dz1dx3n-20sm100_kv2e-4' %14
...
        % 'topo4_H500Lx3k_s0.2dz1dx3n-20sm100_kv1e-4' %40
        % 'topo4_H500Lx3k_s0.4dz1dx3n-20sm100_kv1e-4' %40
        % 'topo4_H500Lx3k_s0.6dz1dx3n-20sm100_kv1e-4' %40
        % 'topo4_H500Lx3k_s0.8dz1dx3n-20sm100_kv1e-4' %32
        % 'topo4_H500Lx3k_s1.0dz1dx3n-20sm100_kv1e-4' %25
        % 'topo4_H500Lx3k_s1.2dz1dx3n-20sm100_kv1e-4'
        % 'topo4_H500Lx3k_s1.4dz1dx3n-20sm100_kv1e-4'
        % 'topo4_H500Lx3k_s1.6dz1dx3n-20sm100_kv1e-4'
        % 'topo4_H500Lx3k_s1.8dz1dx3n-20sm100_kv1e-4'
        ...
        % 'topo4_H500Lx3k_s0.2dz1dx3n-20sm100_kv1e-5' 
        % 'topo4_H500Lx3k_s0.4dz1dx3n-20sm100_kv1e-5'
        % 'topo4_H500Lx3k_s0.6dz1dx3n-20sm100_kv1e-5'
...
        % 'topo4_H500Lx3k_s0.2dz1dx3n-20sm100_kv5e-5'
        % 'topo4_H500Lx3k_s0.4dz1dx3n-20sm100_kv5e-5'
        % 'topo4_H500Lx3k_s0.6dz1dx3n-20sm100_kv5e-5'
        ...
        ...

        % 'topo4_H600Lx10k_s0.0dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s0.1dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s0.2dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s0.3dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s0.4dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s0.5dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s0.6dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s0.7dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s0.8dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s0.9dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s1.0dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s1.1dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s1.2dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s1.3dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s1.4dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s1.5dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s1.6dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s1.7dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s1.8dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s1.9dz2dx10n-20sm50_kv1e-4'
        % 'topo4_H600Lx10k_s2.0dz2dx10n-20sm50_kv1e-4'
    ...
    %     'topo0_H600Lx10k_s0.0dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s0.1dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s0.2dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s0.3dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s0.4dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s0.5dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s0.6dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s0.7dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s0.8dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s0.9dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s1.0dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s1.1dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s1.2dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s1.3dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s1.4dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s1.5dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s1.6dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s1.7dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s1.8dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s1.9dz2dx10n-20sm50_kv1e-4'
    %     'topo0_H600Lx10k_s2.0dz2dx10n-20sm50_kv1e-4'
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
% 'topo0_H500_s0.002dz1dx3ln200n-20sm100_kv1e-4'
% 'topo0_H500_s0.002dz1dx3ln200n-20sm100_kv1e-4_smag'
% 'topo0_H500_s0.002dz1dx3ln200n-20sm100_kv2e-4_restoreT'
% 'topo0_H500_s0.002dz1dx3ln200n-20sm100_kv2e-4_smag'
                ...
% 'topo0_H500_s0dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0001dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0002dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0003dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0004dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0005dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0006dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0007dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0008dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0009dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.001dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0011dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0012dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0013dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0014dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0015dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0016dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0017dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0018dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.0019dz1dx3ln200n-20sm100_kv2e-4'
% 'topo0_H500_s0.002dz1dx3ln200n-20sm100_kv2e-4'
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
        ...
        % 'linear_topo0_H500_s0.001dz1dx5ln200n-10s_kv2e-4'
        % 'linear_topo0_H500_s0.0015dz1dx5ln200n-10s_kv2e-4'
        % 'linear_topo0_H500_s0.002dz1dx5ln200n-10s_kv2e-4'
        % 'linear_topo0_H500_s0.0025dz1dx5ln200n-10s_kv2e-4'
        % 'norestore_linear_topo0_H500_s0.001dz1dx5ln200n-10s_kv2e-4'
        % 'norestore_linear_topo0_H500_s0.0015dz1dx5ln200n-10s_kv2e-4'
        % 'norestore_linear_topo0_H500_s0.002dz1dx5ln200n-10s_kv2e-4'
        % 'norestore_linear_topo0_H500_s0.0025dz1dx5ln200n-10s_kv2e-4'
        % ...
        % 'tanh_topo0_H500_s0.001dz1dx5ln200n-10s_kv2e-4'
        % 'tanh_topo0_H500_s0.0015dz1dx5ln200n-10s_kv2e-4'
        % 'tanh_topo0_H500_s0.002dz1dx5ln200n-10s_kv2e-4'
        % 'tanh_topo0_H500_s0.0025dz1dx5ln200n-10s_kv2e-4'
        % 'norestore_tanh_topo0_H500_s0.001dz1dx5ln200n-10s_kv2e-4'
        % 'norestore_tanh_topo0_H500_s0.0015dz1dx5ln200n-10s_kv2e-4'
        % 'norestore_tanh_topo0_H500_s0.002dz1dx5ln200n-10s_kv2e-4'
        % 'norestore_tanh_topo0_H500_s0.0025dz1dx5ln200n-10s_kv2e-4'
        % ...
        % 'tanh2_topo0_H500_s0.001dz1dx5ln200n-10s_kv2e-4'
        % 'tanh2_topo0_H500_s0.0015dz1dx5ln200n-10s_kv2e-4'
        % 'tanh2_topo0_H500_s0.002dz1dx5ln200n-10s_kv2e-4'
        % 'tanh2_topo0_H500_s0.0025dz1dx5ln200n-10s_kv2e-4'
        % 'norestore_tanh2_topo0_H500_s0.001dz1dx5ln200n-10s_kv2e-4'
        % 'norestore_tanh2_topo0_H500_s0.0015dz1dx5ln200n-10s_kv2e-4'
        % 'norestore_tanh2_topo0_H500_s0.002dz1dx5ln200n-10s_kv2e-4'
        % 'norestore_tanh2_topo0_H500_s0.0025dz1dx5ln200n-10s_kv2e-4'
        ...
        % 'Lx5kmRestoreT_topo0_H500_s0.0007dz1dx10ln200n-14s_kv2e-4'
        ...
        % 'Lx10kmRestoreT_topo0_H500_s0dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0001dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0002dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0003dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0004dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0005dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0006dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0007dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0008dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0009dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.001dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0011dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0012dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0013dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0014dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0015dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0016dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0017dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0018dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.0019dz1dx10ln200n-20_kv2e-4'
        % 'Lx10kmRestoreT_topo0_H500_s0.002dz1dx10ln200n-20_kv2e-4'
        ...
        % 'topo0_tanhH800_s0dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0001dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0002dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0003dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0004dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0005dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0006dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0007dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0008dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0009dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.001dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0011dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0012dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0013dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0014dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0015dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0016dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0017dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0018dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.0019dz1dx3ln200n-20_kv2e-4'
        % 'topo0_tanhH800_s0.002dz1dx3ln200n-20_kv2e-4'




%         'test_lores_topo0_tanhH800_s0.0003dz3dx20ln200n-20_kv2e-4'
%         'test_lores_topo0_tanhH800_s0.0009dz3dx20ln200n-20_kv2e-4'
%         'test_lores_topo0_tanhH800_s0.0018dz3dx20ln200n-20_kv2e-4'

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

    % ...
% 'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_nosmag_h1e-4v2e-4'
% 'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_nosmag_hv2e-4'
% 'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_nosmag_h4e-4v2e-4'
% % 'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_nosmag_h2e-4v1e-4'
% % 'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_nosmag_h6e-4v3e-4'
% % 'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_nosmag_h1e-3v1e-4'
% ... 
% 'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_smag_h1e-4v2e-4'
% 'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_smag_h1e-4v5e-4'
% % 'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_smag_h1e-4v1e-3'
% 'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_smag_h5e-4v2e-4'
% ...
% 'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_hv1e-5_3DSmag1e-2'
% 'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_hv1e-5_3DSmag1e-3'
% 'topo0_H900_smo90m_s0.0006_dz3dx20_no1e-7ln200_hv1e-5_3DSmag1e-4' %%% Too small     
% ...
% 'topo0_H900_smo90m_s0.0015_dz3dx20_no1e-7ln200_hv1e-5_3DSmag1e-2' %%% Significantly reduce the growth rate
% 'topo0_H900_smo90m_s0.0015_dz3dx20_no1e-7ln200_hv1e-5_3DSmag1e-3'
% 'topo0_H900_smo90m_s0.0015_dz3dx20_no1e-7ln200_hv1e-5_3DSmag1e-4'
% ...
% 'topo0_H900_smo90m_s0.0012_dz3dx20_no1e-7ln200_NOsmag_h1e-4v2e-4'
% 'topo0_H900_smo90m_s0.0012_dz3dx20_no1e-7ln200_NOsmag_h1e-4v5e-4'
% 'topo0_H900_smo90m_s0.0012_dz3dx20_no1e-7ln200_NOsmag_h2e-4v2e-4'
% 'topo0_H900_smo90m_s0.0012_dz3dx20_no1e-7ln200_smag_h1e-4v2e-4'
% 'topo0_H900_smo90m_s0.0012_dz3dx20_no1e-7ln200_smag_h1e-4v5e-4'
% 'topo0_H900_smo90m_s0.0012_dz3dx20_no1e-7ln200_smag_h2e-4v2e-4'
% ...
% 'topo0_H900_smo90m_s0.0015_dz3dx20_no1e-7ln200_NOsmag_h1e-4v2e-4'
% 'topo0_H900_smo90m_s0.0015_dz3dx20_no1e-7ln200_NOsmag_h1e-4v5e-4'
% 'topo0_H900_smo90m_s0.0015_dz3dx20_no1e-7ln200_NOsmag_h2e-4v2e-4'
% 'topo0_H900_smo90m_s0.0015_dz3dx20_no1e-7ln200_smag_h1e-4v2e-4'
% 'topo0_H900_smo90m_s0.0015_dz3dx20_no1e-7ln200_smag_h1e-4v5e-4'
% 'topo0_H900_smo90m_s0.0015_dz3dx20_no1e-7ln200_smag_h2e-4v2e-4'


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