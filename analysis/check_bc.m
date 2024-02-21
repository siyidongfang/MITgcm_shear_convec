%%%
%%% check_bc.m
%%%
%%% Check the bottom boundary conditions

ne = 1;
load_all

No = nDumps-1;

for o=No-10

    nIter = dumpIters(o);
    time_h = nIter.*deltaT./3600;
    time_tidal(o) = time_h/12;

    tt = squeeze(rdmds([exppath,'/results/THETA'],nIter));
    uu = squeeze(rdmds([exppath,'/results/UVEL'],nIter));
    % vv = squeeze(rdmds([exppath,'/results/VVEL'],nIter));
    ww = squeeze(rdmds([exppath,'/results/WVEL'],nIter));

end

figure(1)
pcolor(tt);shading flat;colorbar;colormap(redblue);clim([-0.1 0.1])

