

    expdir ='../../MITgcm_shear_convec/exps_topo4_kv2e-5/'
    % expdir ='../../MITgcm_shear_convec/exps_flat_noCori/'
    % expdir ='../../MITgcm_shear_convec/exps_topo4_noCori/'
    % expdir ='../../MITgcm_shear_convec/exps_topo4_tanh/'
    % expdir ='../../MITgcm_shear_convec/exps_topo4_tanh_zeroCenter/'

    % NTtide_all = [100 100 100 72 64 56 49 45];
    % Nstart_all = [50  83  20  30 10 20 30 30];
    % Nend_all   = [100 100 30  40 20 30 49 45];

    EXPNAME = {
        'topo4_H500Lx3k_s0.8dz1dx3n-20sm100_kv1e-6' %32
        'topo4_H500Lx3k_s0.8dz1dx3n-20sm100_kv2e-6' %32
        'topo4_H500Lx3k_s0.8dz1dx3n-20sm100_kv5e-6' %32
        'topo4_H500Lx3k_s0.8dz1dx3n-20sm100_kv7.5e-6' %32
        'topo4_H500Lx3k_s0.8dz1dx2.1873n-20sm100_kv1e-5' %23
        'topo4_H500Lx3k_s0.8dz1dx2.1873n-20sm100_kv8e-5' %23
        'topo4_H500Lx3k_s0.8dz1dx6.562n-20sm100_kv8e-5'  %68
        ...
        'topo4_H500Lx3k_s0.5dz1dx3n-20sm100_kv2e-5' %49
        'topo4_H500Lx3k_s0.8dz1dx3n-20sm100_kv1e-5' %32
        'topo4_H500Lx3k_s0.8dz1dx3n-20sm100_kv2e-5' %31
        'topo4_H500Lx3k_s0.8dz1dx3n-20sm100_kv4e-5' %31
        'topo4_H500Lx3k_s1.0dz1dx3n-20sm100_kv2e-5' %25
        'topo4_H500Lx3k_s1.4dz1dx3n-20sm100_kv2e-5' %18
        ...
        % 'topo4_H500Lx3k_s0.3dz1dx3n-20sm100_kv6.5e-5' %82
        % 'topo4_H500Lx3k_s0.5dz1dx3n-20sm100_kv6.5e-5' %48
        % 'topo4_H500Lx3k_s0.8dz1dx3n-20sm100_kv6.5e-5' %31
        % 'topo4_H500Lx3k_s1.0dz1dx3n-20sm100_kv6.5e-5' %25
        % 'topo4_H500Lx3k_s1.2dz1dx3n-20sm100_kv6.5e-5' %21
        % 'topo4_H500Lx3k_s1.4dz1dx3n-20sm100_kv6.5e-5' %18
        ...
        % 'topo4_H500Lx3k_s0.1dz1dx3n-20sm100_kv8e-5' %100
        % 'topo4_H500Lx3k_s0.3dz1dx3n-20sm100_kv8e-5' %83 74
        % 'topo4_H500Lx3k_s0.5dz1dx3n-20sm100_kv8e-5' %151
        % 'topo4_H500Lx3k_s0.7dz1dx3n-20sm100_kv8e-5' %108 
        % 'topo4_H500Lx3k_s0.8dz1dx3n-20sm100_kv8e-5' %96 29
        % 'topo4_H500Lx3k_s0.9dz1dx3n-20sm100_kv8e-5' %84 25
        % 'topo4_H500Lx3k_s1.0dz1dx3n-20sm100_kv8e-5' %75
        % 'topo4_H500Lx3k_s1.1dz1dx3n-20sm100_kv8e-5' %46
        % 'topo4_H500Lx3k_s1.2dz1dx3n-20sm100_kv8e-5' %61
        % 'topo4_H500Lx3k_s1.3dz1dx3n-20sm100_kv8e-5' %19
        % 'topo4_H500Lx3k_s1.4dz1dx3n-20sm100_kv8e-5' %36
        % 'topo4_H500Lx3k_s1.5dz1dx3n-20sm100_kv8e-5' %17
        % 'topo4_H500Lx3k_s1.6dz1dx3n-20sm100_kv8e-5' %16
        % 'topo4_H500Lx3k_s1.7dz1dx3n-20sm100_kv8e-5' %15
        % 'topo4_H500Lx3k_s1.8dz1dx3n-20sm100_kv8e-5' %14
        % 'topo4_H500Lx3k_s1.9dz1dx3n-20sm100_kv8e-5' %13
        % 'topo4_H500Lx3k_s2.0dz1dx3n-20sm100_kv8e-5' %13
        % 'topo4_H500Lx3k_s2.1dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H500Lx3k_s2.2dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H500Lx3k_s2.3dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H500Lx3k_s2.4dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H500Lx3k_s2.5dz1dx3n-20sm100_kv8e-5' %
        ...
        % 'topo0_H500Lx3k_s0.1dz1dx3n-20sm100_kv8e-5' %100
        % 'topo0_H500Lx3k_s0.3dz1dx3n-20sm100_kv8e-5' %83
        % 'topo0_H500Lx3k_s0.5dz1dx3n-20sm100_kv8e-5' %50
        % 'topo0_H500Lx3k_s0.7dz1dx3n-20sm100_kv8e-5' %72
        % 'topo0_H500Lx3k_s0.8dz1dx3n-20sm100_kv8e-5' %64
        % 'topo0_H500Lx3k_s0.9dz1dx3n-20sm100_kv8e-5' %56
        % 'topo0_H500Lx3k_s1.0dz1dx3n-20sm100_kv8e-5' %50
        % 'topo0_H500Lx3k_s1.1dz1dx3n-20sm100_kv8e-5' %45
        % 'topo0_H500Lx3k_s1.2dz1dx3n-20sm100_kv8e-5' %21*2
        % 'topo0_H500Lx3k_s1.3dz1dx3n-20sm100_kv8e-5' %20*2
        % 'topo0_H500Lx3k_s1.4dz1dx3n-20sm100_kv8e-5' %18
        % 'topo0_H500Lx3k_s1.5dz1dx3n-20sm100_kv8e-5' %34
        % 'topo0_H500Lx3k_s1.6dz1dx3n-20sm100_kv8e-5' %16
        % 'topo0_H500Lx3k_s1.7dz1dx3n-20sm100_kv8e-5' %15
        % 'topo0_H500Lx3k_s1.8dz1dx3n-20sm100_kv8e-5' %14
        % 'topo0_H500Lx3k_s1.9dz1dx3n-20sm100_kv8e-5' %13
        % 'topo0_H500Lx3k_s2.0dz1dx3n-20sm100_kv8e-5' %13
        % 'topo0_H500Lx3k_s2.1dz1dx3n-20sm100_kv8e-5' %
        % 'topo0_H500Lx3k_s2.2dz1dx3n-20sm100_kv8e-5' %
        % 'topo0_H500Lx3k_s2.3dz1dx3n-20sm100_kv8e-5' %
        % 'topo0_H500Lx3k_s2.4dz1dx3n-20sm100_kv8e-5' %
        % 'topo0_H500Lx3k_s2.5dz1dx3n-20sm100_kv8e-5' %
        ...
        % 'topo4_H800Lx3k_s0.1dz1dx3n-20sm100_kv8e-5' %64
        % 'topo4_H800Lx3k_s0.3dz1dx3n-20sm100_kv8e-5' %35
        % 'topo4_H800Lx3k_s0.5dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s0.7dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s0.8dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s0.9dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s1.0dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s1.1dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s1.2dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s1.3dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s1.4dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s1.5dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s1.6dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s1.7dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s1.8dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s1.9dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s2.0dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s2.1dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s2.2dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s2.3dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s2.4dz1dx3n-20sm100_kv8e-5' %
        % 'topo4_H800Lx3k_s2.5dz1dx3n-20sm100_kv8e-5' %
       };
    nEXP = length(EXPNAME);



