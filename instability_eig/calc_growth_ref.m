
clear;

constants_ref;
% Ntopo = 1; % flat
Ntopo = 2; % topo4

grow_eig = NaN.*zeros(1,Ns);
max_kidx = NaN.*zeros(1,Ns);
max_midx = NaN.*zeros(1,Ns);
grow_eig_km = NaN.*zeros(Nk,Ns);

for ntopo = Ntopo
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
            grow_eig_km(:,ns) = grow;
            [grow_eig(ns),I] =max(grow,[],'all','omitnan');
            [max_kidx(ns), max_midx(ns)] = ind2sub(size(grow), I);
        end
    end
end

save('products/grow_topo4_eig_km.mat');


function grow =load_func(fname)
         load(fname,'grow');
end