
    % expdir ='../../MITgcm_shear_convec/exps_kv5e-6/'
    expdir ='../../MITgcm_shear_convec/exps_kv1e-5/'

    EXPNAME = {
        'topo4_H500Lx3k_s0.100dz1dx3sm100' %100
        'topo4_H500Lx3k_s0.400dz1dx3sm100' %62
        'topo4_H500Lx3k_s0.600dz1dx3sm100' %42
        'topo4_H500Lx3k_s0.800dz1dx3sm100' %32
        'topo4_H500Lx3k_s1.000dz1dx3sm100' %25.5
        'topo4_H500Lx3k_s1.200dz1dx3sm100' %21.5
        'topo4_H500Lx3k_s1.400dz1dx3sm100' %18.5
        'topo4_H500Lx3k_s1.500dz1dx3sm100' %17.25
        'topo4_H500Lx3k_s1.600dz1dx3sm100' %15.75
        'topo4_H500Lx3k_s1.700dz1dx3sm100' %15.25
        'topo4_H500Lx3k_s1.800dz1dx3sm100' %14.5
        'topo4_H500Lx3k_s1.900dz1dx3sm100' %13.5
        'topo4_H500Lx3k_s2.000dz1dx3sm100' %13
        'topo4_H500Lx3k_s2.070dz1dx3sm100' %12.5
        ...
        % 'tanh_topo4_H800Lx3k_s0.400dz1dx3sm100' %38
        % 'tanh_topo4_H800Lx3k_s1.000dz1dx3sm100' %15.5 -- resubmit
        % 'tanh_topo4_H800Lx3k_s1.500dz1dx3sm100' %10.5
        % 'tanh_topo4_H800Lx3k_s1.700dz1dx3sm100' %9
        % 'tanh_topo4_H800Lx3k_s1.800dz1dx3sm100' %8.75
        % 'tanh_topo4_H800Lx3k_s1.900dz1dx3sm100' %8.25
        % 'tanh_topo4_H800Lx3k_s2.000dz1dx3sm100' %7.75
        % 'tanh_topo4_H800Lx3k_s2.070dz1dx3sm100' %7.5
        % ...
        % 'zeroCenter_tanh_topo4_H800Lx3k_s0.400dz1dx3sm100' %?? No grow
        % 'zeroCenter_tanh_topo4_H800Lx3k_s1.000dz1dx3sm100' %?? nonlinear rapid grow
        % 'zeroCenter_tanh_topo4_H800Lx3k_s1.500dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H800Lx3k_s1.700dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H800Lx3k_s1.800dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H800Lx3k_s1.900dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H800Lx3k_s2.000dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H800Lx3k_s2.070dz1dx3sm100'
        ...
        ...
        'zeroCenter_flat_H800Lx3k_s0.400dz1dx3sm100' %38
        'zeroCenter_flat_H800Lx3k_s1.000dz1dx3sm100' %15.5
        'zeroCenter_flat_H800Lx3k_s1.400dz1dx3sm100' %11.25
        'zeroCenter_flat_H800Lx3k_s1.600dz1dx3sm100' %9.75
        'zeroCenter_flat_H800Lx3k_s1.800dz1dx3sm100' %8.75
        'zeroCenter_flat_H800Lx3k_s2.000dz1dx3sm100' %7.75
        'zeroCenter_flat_H800Lx3k_s2.200dz1dx3sm100' %7
        'zeroCenter_flat_H800Lx3k_s2.400dz1dx3sm100' %6.5
        'zeroCenter_flat_H800Lx3k_s2.600dz1dx3sm100' %6
        'zeroCenter_flat_H800Lx3k_s2.800dz1dx3sm100' %5.5
        ...
        % 'tanh_topo0_H800Lx3k_s0.100dz1dx3sm100' %100
        % 'tanh_topo0_H800Lx3k_s0.400dz1dx3sm100' %39
        % 'tanh_topo0_H800Lx3k_s0.600dz1dx3sm100' %26
        % 'tanh_topo0_H800Lx3k_s0.800dz1dx3sm100' %19
        % 'tanh_topo0_H800Lx3k_s1.000dz1dx3sm100' %15 
        % 'tanh_topo0_H800Lx3k_s1.200dz1dx3sm100' %13 --submitted
        % 'tanh_topo0_H800Lx3k_s1.400dz1dx3sm100' %11 --submitted
        % 'tanh_topo0_H800Lx3k_s1.600dz1dx3sm100' %9.5
        % 'tanh_topo0_H800Lx3k_s1.800dz1dx3sm100' %16.75
        % 'tanh_topo0_H800Lx3k_s2.000dz1dx3sm100' %8
        % 'tanh_topo0_H800Lx3k_s2.200dz1dx3sm100' %7.25
        % 'tanh_topo0_H800Lx3k_s2.400dz1dx3sm100' %6.5
        % 'tanh_topo0_H800Lx3k_s2.600dz1dx3sm100' %12
        % 'tanh_topo0_H800Lx3k_s2.800dz1dx3sm100' %5.75
        % 'flat_H500Lx3k_s0.10dz1dx3sm100' %100
        % 'flat_H500Lx3k_s0.40dz1dx3sm100' %63
        % 'flat_H500Lx3k_s0.60dz1dx3sm100' %43
        % 'flat_H500Lx3k_s0.80dz1dx3sm100' %32
        % 'flat_H500Lx3k_s1.00dz1dx3sm100' %26
        % 'flat_H500Lx3k_s1.10dz1dx3sm100' %23
        % 'flat_H500Lx3k_s1.20dz1dx3sm100' %21
        % 'flat_H500Lx3k_s1.30dz1dx3sm100' %20
        % 'flat_H500Lx3k_s1.40dz1dx3sm100' %19
        % 'flat_H500Lx3k_s1.455dz1dx3sm100'%18
        % 'flat_H500Lx3k_s1.50dz1dx3sm100' %17
        % 'flat_H500Lx3k_s1.60dz1dx3sm100'
        % 'flat_H500Lx3k_s1.70dz1dx3sm100'
        % 'flat_H500Lx3k_s1.80dz1dx3sm100'
        % 'flat_H500Lx3k_s1.90dz1dx3sm100'
        % 'flat_H500Lx3k_s2.00dz1dx3sm100'
        % 'flat_H500Lx3k_s2.10dz1dx3sm100'
        % 'flat_H500Lx3k_s2.20dz1dx3sm100'
        % 'flat_H500Lx3k_s2.400dz1dx3sm100' %11
        % 'flat_H500Lx3k_s2.600dz1dx3sm100' %10
        % 'flat_H500Lx3k_s2.800dz1dx3sm100'
        ...
        % 'topo4_H500Lx3k_s0.10dz1dx3sm100'
        % 'topo4_H500Lx3k_s0.40dz1dx3sm100'
        % 'topo4_H500Lx3k_s0.60dz1dx3sm100'
        % 'topo4_H500Lx3k_s0.80dz1dx3sm100' %64
        % 'topo4_H500Lx3k_s1.00dz1dx3sm100' %52
        % 'topo4_H500Lx3k_s1.10dz1dx3sm100'
        % 'topo4_H500Lx3k_s1.20dz1dx3sm100' % submitted twice
        % 'topo4_H500Lx3k_s1.30dz1dx3sm100'
        % 'topo4_H500Lx3k_s1.335dz1dx3sm100'
        % 'topo4_H500Lx3k_s1.40dz1dx3sm100' % submitted twice
        % 'topo4_H500Lx3k_s1.50dz1dx3sm100'
        % 'topo4_H500Lx3k_s1.60dz1dx3sm100'
        % 'topo4_H500Lx3k_s1.70dz1dx3sm100'
        % 'topo4_H500Lx3k_s1.80dz1dx3sm100'
        % 'topo4_H500Lx3k_s1.90dz1dx3sm100'
        % 'topo4_H500Lx3k_s2.00dz1dx3sm100'
        % 'topo4_H500Lx3k_s2.08dz1dx3sm100'
        % ...
        % 'tanh_topo4_H500Lx3k_s0.10dz1dx3sm100' %100
        % 'tanh_topo4_H500Lx3k_s0.40dz1dx3sm100' %38
        % 'tanh_topo4_H500Lx3k_s0.60dz1dx3sm100' %25 - submitted
        % 'tanh_topo4_H500Lx3k_s0.80dz1dx3sm100' %19 - submitted
        % 'tanh_topo4_H500Lx3k_s1.00dz1dx3sm100' %15 - submitted
        % 'tanh_topo4_H500Lx3k_s1.10dz1dx3sm100' % - submitted
        % 'tanh_topo4_H500Lx3k_s1.20dz1dx3sm100'
        % 'tanh_topo4_H500Lx3k_s1.30dz1dx3sm100'
        % 'tanh_topo4_H500Lx3k_s1.335dz1dx3sm100'
        % 'tanh_topo4_H500Lx3k_s1.40dz1dx3sm100'
        % 'tanh_topo4_H500Lx3k_s1.50dz1dx3sm100'
        % 'tanh_topo4_H500Lx3k_s1.60dz1dx3sm100'
        % 'tanh_topo4_H500Lx3k_s1.70dz1dx3sm100'
        % 'tanh_topo4_H500Lx3k_s1.80dz1dx3sm100'
        % 'tanh_topo4_H500Lx3k_s1.90dz1dx3sm100'
        % 'tanh_topo4_H500Lx3k_s2.00dz1dx3sm100'
        % 'tanh_topo4_H500Lx3k_s2.08dz1dx3sm100'
        ...
        % 'zeroCenter_tanh_topo4_H500Lx3k_s0.10dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s0.40dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s0.60dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s0.80dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s1.00dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s1.10dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s1.20dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s1.30dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s1.335dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s1.40dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s1.50dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s1.60dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s1.70dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s1.80dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s1.90dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s2.00dz1dx3sm100'
        % 'zeroCenter_tanh_topo4_H500Lx3k_s2.08dz1dx3sm100'
       };
    nEXP = length(EXPNAME);



