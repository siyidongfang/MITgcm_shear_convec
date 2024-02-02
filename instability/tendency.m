    


    zeta0 = @(z) interp1(linspace(min(zspan),max(zspan),numel(z0)), z0, z);    
    %%% Form Initial Guess
    xmesh = linspace(0,1,Nr+1);
    if(sum(z0~=0)>0)
    solinit = bvpinit(xmesh, [0.1 0]);
    else
    solinit = bvpinit(xmesh, [0 0]);
    end
    %%% Solve the boundary value problem using the bvp4c solver.
    sol1 = bvp4c(@(z,y)bvpfun(z,y,kx,zeta0), @bcfun, solinit);
    psi0 = sol1.y(1,:);
    zz0 = sol1.x;
    % p0 = interp1(zz0,psi0,zz_wgrid);
    p0 = psi0;


    %%% Boundary condition:  
    % b0(1) = b0(2); 
    % b0(Nr) = b0(Nr-1);

    U = cos(t0)*Atide/U0;
    U_wgrid = cos(t0)*Atide_wgrid/U0;
    if(U0==0)
        U=zeros(1,Nr);
        U_wgrid = zeros(1,Nr+1);
    end

    for m = 2:Nr-1
        d2bdz2(m) = (b0(m-1)-2*b0(m)+b0(m+1))/dz^2;
    end

    dUdz = (U_wgrid(2:Nr+1)-U_wgrid(1:Nr))/dz;
    d2bdz2(1) = 0; d2bdz2(Nr) = 0;

    p0_ugrid = 0.5*(p0(1:Nr)+p0(2:Nr+1));
    dpsidz = (p0(2:Nr+1)-p0(1:Nr))/dz;

    bq1(o,:) = -1i*kx*C1*U.*b0;
    bq2(o,:) = -1i*kx*cotd(topo).*p0_ugrid;
    bq3(o,:) = dpsidz;
    bq4(o,:) = +1i*kx*C1*dUdz.*tan(t0).*p0_ugrid;
    bq5(o,:) = +C4*d2bdz2-C4*kx^2.*b0;
    if(NOdiffusion)
        bq5(o,:) = zeros(1,Nr); %%% ignore diffusion
    end

    dbdt(o,:) = bq1(o,:) + bq2(o,:) + bq3(o,:) ...
              + bq4(o,:) + bq5(o,:);


    dbdz = zeros(1,Nr+1);
    dbdz(2:Nr) = (b0(2:Nr)-b0(1:Nr-1))/dz;
    dbdz(1) = 0; dbdz(Nr+1) = 0;

    b0_wgrid = zeros(1,Nr+1);
    b0_wgrid(2:Nr) = 0.5*(b0(2:Nr)+b0(1:Nr-1));
    b0_wgrid(1) = b0_wgrid(2); b0_wgrid(Nr+1) = b0_wgrid(Nr);

    for m = 2:Nr
        d2zetadz2(m) = (z0(m-1)-2*z0(m)+z0(m+1))/dz^2;
    end
    zeta(1) = 0; zeta(Nr+1) = 0; %%% Free-slip

    % !!! d2zetadz2(1) = d2zetadz2(2); d2zetadz2(Nr+1) = d2zetadz2(Nr); %????? Is this correct?

    zq1(o,:) = -1i*kx*C1*U_wgrid.*z0;
    zq2(o,:) = +C2^2*(1i*kx*cotd(topo)*b0_wgrid);
    zq3(o,:) = +C2^2*(-dbdz);
    zq4(o,:) = +C3*d2zetadz2-C3*kx^2.*z0;
    if(NOdiffusion)
        zq4(o,:) = zeros(1,Nr+1); %%% ignore dissipation
    end

    dzetadt(o,:) = zq1(o,:) + zq2(o,:) + zq3(o,:) + zq4(o,:);


    %%% Code Equation
    function dydz = bvpfun(z,y,kx,zeta)
        dydz = [y(2); kx^2*y(1)+zeta(z)];
    end
    
    %%% Code Boundary Conditions
    function res = bcfun(ya,yb)
        res = [ya(1)
               yb(1)];
    end


