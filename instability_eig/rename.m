


constants;

for ntopo = 2:Ntopo
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
            oldname = [expdir 'ptide' num2str(Ptide) '_topo' num2str(topo) '_N' num2str(N) '_shear' num2str(shear*1e3,3) '.mat'];
            newname = [expdir 'ptide' num2str(Ptide) '_topo' num2str(topo) '_N' num2str(N*1e3,3) '_shear' num2str(shear*1e3,3) '.mat'];
            movefile(oldname, newname); 
        end
    end
end



