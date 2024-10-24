clear;close all;


load('products/grow_M2_calc.mat','Ri','omega','Ns','kx_all','m0_rw','Nk','topo_all','N_all')

for ntopo = 1:8
    ntopo
for nn = 3

        clear shear_all Ri_shear omega0 grow omega0_sort epsilon_sort grow_sort Ri_sort

        topo=topo_all(ntopo);
        N=N_all(nn);
        shear_all = [0:3*N/(Ns-1):3*N];
        Ri_shear = squeeze(Ri(ntopo,3,:))';
        omega0 = abs(N*kx_all/m0_rw);
        
        grow_all = [];
        epsilon_all = [];
        omega0_all = [];
        Ri_all = [];
        
        
        for ns=1:Ns
            shear = shear_all(ns);
            load(['exps_sens_M2/ptide43200_topo' num2str(topo) '_N' num2str(N*1e3) '_shear' num2str(shear*1e3,'%.4f') '.mat'],'grow');
            grow(grow<=0)=NaN;
            grow = log(grow)/43200*3600;
        
            grow_all = [grow_all grow];
            epsilon = 2*omega0.^3/omega*shear/N;
            epsilon_all=[epsilon_all epsilon];
            omega0_all = [omega0_all omega0];
            Ri_all = [Ri_all Ri_shear(ns).*ones(1,length(grow))];
        end
        
            [omega0_sort,I] = sort(omega0_all);
            epsilon_sort = epsilon_all(I);
            grow_sort = grow_all(I);
            Ri_sort = Ri_all(I);
        
        
        save(['products/methieu_M2_topo' num2str(topo) '_N' num2str(N*1e3) '.mat'],'topo','N','Ri_sort','grow_sort','epsilon_sort','omega0_sort')



end
end