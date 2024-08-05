

    expdir ='../../MITgcm_shear_convec/exps_flat_hires/'
    % expdir ='../../MITgcm_shear_convec/exps_topo4_hires/'

    EXPNAME = {
        'topo4_H500Lx3k_s0.0dz1dx3n-20sm100_kv8e-5' %20
        'topo4_H500Lx3k_s0.1dz1dx3n-20sm100_kv8e-5' %100
        'topo4_H500Lx3k_s0.2dz1dx3n-20sm100_kv8e-5' %100
        'topo4_H500Lx3k_s0.3dz1dx3n-20sm100_kv8e-5' %83
        'topo4_H500Lx3k_s0.4dz1dx3n-20sm100_kv8e-5' %63
        'topo4_H500Lx3k_s0.5dz1dx3n-20sm100_kv8e-5' %49
        'topo4_H500Lx3k_s0.6dz1dx3n-20sm100_kv8e-5' %42
        'topo4_H500Lx3k_s0.7dz1dx3n-20sm100_kv8e-5' %36
        'topo4_H500Lx3k_s0.8dz1dx3n-20sm100_kv8e-5' %32
        'topo4_H500Lx3k_s0.9dz1dx3n-20sm100_kv8e-5' %28
        'topo4_H500Lx3k_s1.0dz1dx3n-20sm100_kv8e-5' %26
        'topo4_H500Lx3k_s1.1dz1dx3n-20sm100_kv8e-5' %23
        'topo4_H500Lx3k_s1.2dz1dx3n-20sm100_kv8e-5' %21
        'topo4_H500Lx3k_s1.3dz1dx3n-20sm100_kv8e-5' %19
        'topo4_H500Lx3k_s1.4dz1dx3n-20sm100_kv8e-5' %18
        'topo4_H500Lx3k_s1.5dz1dx3n-20sm100_kv8e-5' %17
        'topo4_H500Lx3k_s1.6dz1dx3n-20sm100_kv8e-5' %16
        'topo4_H500Lx3k_s1.7dz1dx3n-20sm100_kv8e-5' %15
        'topo4_H500Lx3k_s1.8dz1dx3n-20sm100_kv8e-5' %14
        'topo4_H500Lx3k_s2.0dz1dx3n-20sm100_kv8e-5' %13
        ...
        'topo0_H500Lx3k_s0.0dz1dx3n-20sm100_kv8e-5' %20
        'topo0_H500Lx3k_s0.1dz1dx3n-20sm100_kv8e-5' %100
        'topo0_H500Lx3k_s0.2dz1dx3n-20sm100_kv8e-5' %80
        'topo0_H500Lx3k_s0.3dz1dx3n-20sm100_kv8e-5' %83
        'topo0_H500Lx3k_s0.4dz1dx3n-20sm100_kv8e-5' %63
        'topo0_H500Lx3k_s0.5dz1dx3n-20sm100_kv8e-5' %50
        'topo0_H500Lx3k_s0.6dz1dx3n-20sm100_kv8e-5' %42
        'topo0_H500Lx3k_s0.7dz1dx3n-20sm100_kv8e-5' %36
        'topo0_H500Lx3k_s0.8dz1dx3n-20sm100_kv8e-5' %32
        'topo0_H500Lx3k_s0.9dz1dx3n-20sm100_kv8e-5' %28
        'topo0_H500Lx3k_s1.0dz1dx3n-20sm100_kv8e-5' %26
        'topo0_H500Lx3k_s1.1dz1dx3n-20sm100_kv8e-5' %23
        'topo0_H500Lx3k_s1.2dz1dx3n-20sm100_kv8e-5' %21
        'topo0_H500Lx3k_s1.3dz1dx3n-20sm100_kv8e-5' %20
        'topo0_H500Lx3k_s1.4dz1dx3n-20sm100_kv8e-5' %18
        'topo0_H500Lx3k_s1.5dz1dx3n-20sm100_kv8e-5' %17
        'topo0_H500Lx3k_s1.6dz1dx3n-20sm100_kv8e-5' %16
        'topo0_H500Lx3k_s1.7dz1dx3n-20sm100_kv8e-5' %15
        'topo0_H500Lx3k_s1.8dz1dx3n-20sm100_kv8e-5' %14
        'topo0_H500Lx3k_s1.9dz1dx3n-20sm100_kv8e-5' %13
        'topo0_H500Lx3k_s2.0dz1dx3n-20sm100_kv8e-5' %13
        'topo0_H500Lx3k_s2.3dz1dx3n-20sm100_kv8e-5' %11
        ...
        'topo0_H500Lx3k_s1.0dz1dx3n-20sm100_kv2e-4'        %49
        'topo0_H500Lx2.775k_s1.0dz1dx1.85n-20sm100_kv2e-4' %15
        % 'topo0_H500Lx3.08k_s1.0dz1dx3.08n-20sm100_kv2e-4'  %25
        % 'topo0_H500Lx3.48k_s1.0dz1dx3.48n-20sm100_kv2e-4'  %30
        'topo0_H500Lx3.08k_s1.0dz1dx3.08n-20sm100_kv1e-4' %52
        'topo0_H500Lx3.08k_s1.0dz1dx3.08n-20sm100_kv7.5e-5'  %50
        'topo0_H500Lx3.08k_s1.0dz1dx3.08n-20sm100_kv5e-5' %25
       };
    nEXP = length(EXPNAME);



