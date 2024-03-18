    

    %%%%%%%%%%%% B.C.-6 %%%%%%%%%%%%
    z0(1) = 0; z0(Nr+1) = 0; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    zeta0 = @(z) interp1(linspace(min(zspan),max(zspan),numel(z0)), z0, z);    
    %%% Form Initial Guess
    xmesh = zz_wgrid;
    solinit = bvpinit(xmesh, [0 0]);
    % p0_guess = real(p0);
    % solinit.y(1,:) = p0_guess;
    % solinit.y(2,2:end-1) = (p0_guess(3:end)-p0_guess(1:end-2))/dz/2; %%% Centered difference
    % solinit.y(2,1)=(p0_guess(2)-p0_guess(1))/dz;
    % solinit.y(2,end)=(p0_guess(end)-p0_guess(end-1))/dz;

    %%% Solve the boundary value problem using the bvp4c solver.
    % options = bvpset('RelTol', 1e-3, 'AbsTol', 1e-3, 'NMax', 20000);
    options = bvpset('NMax', 10000);
    sol1 = bvp4c(@(z,y)bvpfun(z,y,kx,zeta0), @bcfun, solinit,options);
    % sol1 = bvp4c(@(z,y)bvpfun(z,y,kx,zeta0), @bcfun, solinit);
    psi0 = sol1.y(1,:);
    zz0 = sol1.x;
    if(length(psi0)>length(p0))
        p0 = interp1(zz0,psi0,zz_wgrid);
        % warning('length(psi0)>length(p0)')
    else
        p0 = psi0;
    end

    %%%%%%%%%%%% B.C.-7 %%%%%%%%%%%%
    %%% Impermeable (w=0) At the ocean bottom
    p0(1) = 0;
    %%% Impermeable (w=0) At the upper boundary 
    p0(Nr+1) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    U = cos(t0)*Atide/U0;
    U_wgrid = cos(t0)*Atide_wgrid/U0;
    if(U0==0)
        U=zeros(1,Nr);
        U_wgrid = zeros(1,Nr+1);
    end

    dUdz = (U_wgrid(2:Nr+1)-U_wgrid(1:Nr))/dz;

    dbdz(o,2:Nr) = (b0(2:Nr)-b0(1:Nr-1))/dz; %%% on w-grid

    %%%%%%%%%%%% B.C.-8 %%%%%%%%%%%%
    %%% No buoyancy flux, Ocean bottom, when diffusivity is not zero
    %%% dbdz+dB/dz+dB0/dz = 0 ==>  dbdz = -(dB/dz+dB0/dz)
    dbdz(o,1) = - (omega/Shear) * (delta/Hshear) * (cosd(topo)/sind(topo)) ...
        -(delta/Hshear)*sin(t0); 
    %%% No buoyancy flux, Upper boundary, when diffusivity is not zero
    dbdz(o,Nr+1) = - (omega/Shear) * (delta/Hshear) * (cosd(topo)/sind(topo)) ...
        -(delta/Hshear)*sin(t0); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    b0_wgrid = zeros(1,Nr+1);
    b0_wgrid(2:Nr) = 0.5*(b0(2:Nr)+b0(1:Nr-1));

    %%%%%%%%%%%% B.C.-9 %%%%%%%%%%%%
    b0_wgrid(1) = 2*b0(1)-b0_wgrid(2);
    b0_wgrid(Nr+1) = 2*b0(Nr)-b0_wgrid(Nr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for m = 2:Nr-1
        d2bdz2(o,m) = (b0(m-1)-2*b0(m)+b0(m+1))/dz^2;
    end

    %%%%%%%%%%%% B.C.-10 %%%%%%%%%%%%
    dbdz1 = 0.5*(dbdz(o,1)+dbdz(o,2));
    dbdzNr = 0.5*(dbdz(o,Nr)+dbdz(o,Nr+1));

    bw2 = b0_wgrid(2); 
    b1 = b0(1);
    zw2 = zz_wgrid(2);
    z1 = zz(1);

    bwNr = b0_wgrid(Nr); 
    bNr = b0(Nr);
    zwNr = zz_wgrid(Nr);
    zNr = zz(Nr);

    % Taylor expansion: bw2 ~= b1 + dbdz1*(zw2-z1) + 0.5*d2bdz2(1)*(zw2-z1)^2
    d2bdz2(o,1) = ( bw2 - b1 - dbdz1*(zw2-z1) ) ./ ( 0.5*(zw2-z1)^2 );
    % Taylor expansion: bwNr ~= bNr + dbdzNr*(zwNr-zNr) + 0.5*d2bdz2(Nr)*(zwNr-zNr)^2
    d2bdz2(o,Nr) = ( bwNr - bNr - dbdzNr*(zwNr-zNr) ) ./ ( 0.5*(zwNr-zNr)^2 );

    %%% Linear extrapolation
    % d2bdz2(o,[1 Nr]) = interp1(zz(2:Nr-1),d2bdz2(o,2:Nr-1),zz([1 Nr]),'linear','extrap');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p0_ugrid = 0.5*(p0(1:Nr)+p0(2:Nr+1));
    dpsidz = (p0(2:Nr+1)-p0(1:Nr))/dz;

    bq1(o,:) = -1i*kx*C1*U.*b0;
    bq2(o,:) = -1i*kx*cotd(topo).*p0_ugrid;
    % bq3(o,:) = dpsidz;
    % bq4(o,:) = +1i*kx*C1*dUdz.*tan(t0).*p0_ugrid;
    bq3(o,:) = 0;
    bq4(o,:) = 0;
    bq5(o,:) = +C4*d2bdz2(o,:)-C4*kx^2.*b0;
    
    dbdt(o,:) = bq1(o,:) + bq2(o,:) + bq3(o,:) ...
              + bq4(o,:) + bq5(o,:);


    for m = 2:Nr
        d2zetadz2(o,m) = (z0(m-1)-2*z0(m)+z0(m+1))/dz^2;
    end

    %%%%%%%%%%%% B.C.-11 %%%%%%%%%%%%
    %%% Linear extrapolation
    d2zetadz2(o,:) = interp1(zz_wgrid(2:Nr),d2zetadz2(o,2:Nr),zz_wgrid,'linear','extrap');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    zq1(o,:) = -1i*kx*C1*U_wgrid.*z0;
    zq2(o,:) = +C2^2*(1i*kx*cotd(topo)*b0_wgrid);
    % zq3(o,:) = +C2^2*(-dbdz(o,:));
    zq3(o,:) = 0;
    zq4(o,:) = +C3*d2zetadz2(o,:)-C3*kx^2.*z0;

    dzetadt(o,:) = zq1(o,:) + zq2(o,:) + zq3(o,:) + zq4(o,:);

    %%% Code equations: see https://www.mathworks.com/help/matlab/math/solve-bvp-with-two-solutions.html
    function dydz = bvpfun(z,y,kx,zeta)
        dydz = [y(2)
            kx^2*y(1)+zeta(z)];
    end
    
    %%%%%%%%%%%% B.C.-12 %%%%%%%%%%%%
    %%% Code boundary conditions for the streamfunction: phi = 0 at z=0 and z=1
    function res = bcfun(ya,yb)
        res = [ya(1)
               yb(1)];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    



