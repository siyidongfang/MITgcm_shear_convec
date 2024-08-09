
clear;
constants;

Ns = 31;
grow_floquet = NaN.*zeros(Ntopo,Nn,Ns);
max_kidx = NaN.*zeros(Ntopo,Nn,Ns);
max_midx = NaN.*zeros(Ntopo,Nn,Ns);

for ntopo = 1:Ntopo
    topo = topo_all(ntopo)
    cs = cosd(topo);
    ss = sind(topo);
    for nn = 1:Nn
        N = N_all(nn)
        clear shear_all;
        shear_all = [0:3*N/30:3*N]; 
        Ns = length(shear_all);
        for ns =1:Ns
            shear = shear_all(ns);
            fname = [expdir 'ptide' num2str(Ptide) '_topo' num2str(topo) '_N' num2str(N*1e3,3) '_shear' num2str(shear*1e3,3) '.mat'];
            load(fname,'grow'); 
            % Find the most unstable mode and the maximum growth rate
            [grow_floquet(ntopo,nn,ns),I] =max(grow,[],'all');
            [max_kidx(ntopo,nn,ns), max_midx(ntopo,nn,ns)] = ind2sub(size(grow), I);
        end
    end
end

save([expdir 'grow_M2.mat']);

