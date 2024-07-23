%%%
%%% calc_nonlinear.m
%%% 
%%% Calculate the momentum advection: 
%%% linear advection (w'U_z+Uu'_x)
%%% non-linear advection: (u'u'_x+w'u'_z), (u'v'_x+w'v'_z), (u'w'_x+w'w'_z)

clear;

%%% Load the data
ne=1;
load_all;

zidx = 300:500;

% No=nDumps;
No=48;

time_tidal = NaN*zeros(1,No);
dx = delX(1);
dz = delR(end);

for o=No

    nIter = dumpIters(o);
    time_tidal(o) = nIter.*deltaT./3600/12;

    uu = squeeze(rdmds([exppath,'/results/UVEL'],nIter));
    ww = squeeze(rdmds([exppath,'/results/WVEL'],nIter));
    vv = squeeze(rdmds([exppath,'/results/VVEL'],nIter));

    U = mean(uu(:,zidx));
    V = mean(vv(:,zidx));
    u=uu(:,zidx)-U;
    v=vv(:,zidx)-V;
    w=ww(:,zidx);

    %%% Grid conversion 
    Uz = % on u-grid, w-lev
    % on u-grid, t-lev
    w_grid = 0.5*(w+w([2:Nx 1],:)); % on u-grid, w-lev
    w_grid(:,1:end-1) = 0.5*(w_grid(:,2:end)+w_grid(:,1:end-1)); % on u-grid, t-lev
    w_grid(:,end)=0.5(0+w_grid(:,end));
    

    ux =

    du_dx = (uu([2:Nx 1],:,:) - uu(1:Nx,:,:)) ./ DX;
    du_dz(:,:,2:Nr) = - diff(uu,1,3) ./ DZC; %%% N.B. Does not account for partial cells
    dv_dx = (vv(1:Nx,:,:) - vv([Nx 1:Nx-1],:,:)) ./ DX;
    dv_dz = zeros(Nx,Ny,Nr);
    dv_dz(:,:,2:Nr) = - diff(vv,1,3) ./ DZC; %%% N.B. Does not account for partial cells


    wUz = Uz.*w_grid;
    Uux = 

end

















