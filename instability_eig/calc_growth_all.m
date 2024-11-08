
clear;

constants_sens_K1;

grow_eig = NaN.*zeros(Ntopo,Nn,Ns);
max_kidx = NaN.*zeros(Ntopo,Nn,Ns);
max_midx = NaN.*zeros(Ntopo,Nn,Ns);

for ntopo = 1:Ntopo
    topo = topo_all(ntopo)
    cs = cosd(topo);
    ss = sind(topo);
    for nn = 1:Nn
        N = N_all(nn);
        clear shear_all;
        shear_all = [0:3*N/(Ns-1):3*N]; 
        Ns = length(shear_all);
        parfor ns =1:Ns
            shear = shear_all(ns);
            fname = [expdir 'ptide' num2str(Ptide) '_topo' num2str(topo) '_N' num2str(N*1e3) '_shear' num2str(shear*1e3,'%.4f') '.mat'];
            grow = load_func(fname); 
            grow(grow<=0)=NaN;
            grow = log(grow)/43200*3600;
            % Find the most unstable mode and the maximum growth rate
            [grow_eig(ntopo,nn,ns),I] =max(grow,[],'all','omitnan');
            [max_kidx(ntopo,nn,ns), max_midx(ntopo,nn,ns)] = ind2sub(size(grow), I);
        end
    end
end

save('products/grow_K1_new.mat');


function grow =load_func(fname)
         load(fname,'grow');
end