%%%
%%% BC.m
%%%
%%% Boundary condition

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% At the ocean bottom %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Impermeable (w=0)
psi(o,1) = 0;

%%% No-stress (free-slip) b.c. (du/dz = 0)
zeta(o,1) = delta/Hshear*cos(t0);

%%% No buoyancy flux
dbdz(1) = - (omega/Shear) * (delta/Hshear) * (cosd(topo)/sind(topo)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% At the upper boundary %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Impermeable (w=0)
psi(o,Nr+1) = 0;

%%% No-stress (free-slip) b.c. (du/dz = 0)
zeta(o,Nr+1) = 0;

%%% No buoyancy flux
dbdz(Nr+1) = 1/cosd(topo) * (...
    - (omega/Shear) * (delta/Hshear) / sind(topo) ...
    - 1i * kx * b0(Nr+1) * sind(topo) ...
    );







