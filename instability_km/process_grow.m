clear;close all
% constants;
expdir = 'exps_new/topo4_nu0/'
load([expdir 'output.mat']);

Ns=length(shear_all);
grow_rw = zeros(Ns,Nrw);

for s = 1:Ns
    s
    shear = shear_all(s);

    parfor j=1:Nrw
        fname = [[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '_rw' num2str(j) '.mat'];
        grow = load_func(fname);
        if(length(grow)>1)
            grow_rw(s,j)=grow(j);
        else
            grow_rw(s,j)=grow;
        end
    end
end

save('exps_new/topo4_nu0/grow_rw.mat')


function grow = load_func(file)
    load( file );
end