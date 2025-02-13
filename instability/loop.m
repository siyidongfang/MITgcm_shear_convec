
function [p0,dzetadt,dbdt,bq1,bq2,bq3,bq4,bq5,zq1,zq2,zq3,zq4] ...
            = loop(o,Nr,t0,zspan,zz_wgrid,z0,dz,p0,b0,Atide,Atide_wgrid,dzetadt,dbdt,bq1,bq2,bq3,bq4,bq5,zq1,zq2,zq3,zq4,...
                   kx,omega,topo,nu,kappa,N,dAdz,noBQ2,noBQ3,noBQ4,noZQ2,noZQ3,hydrostatic,useTanhShear)

    %%%%%%%%%%%% B.C.-6 %%%%%%%%%%%%
    z0(1) = 0; z0(Nr+1) = 0; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Solve the boundary value problem using the bvp4c solver.
    zeta0 = @(z) interp1(linspace(min(zspan),max(zspan),numel(z0)), z0, z);    
    %%%%% Form Initial Guess -- old
    options = bvpset('NMax', 15000);
    xmesh = zz_wgrid;
    solinit = bvpinit(xmesh, [0 0]);
    %%%%% Form Initial Guess -- new
    % options = bvpset('RelTol', 1e-2, 'AbsTol', 100, 'NMax', 20000);
    % xmesh = zz_wgrid;
    % solinit = bvpinit(xmesh,@(x)guess(x,max(zspan)));
    %%%%% Form Initial Guess -- new
    % options = bvpset('RelTol', 1e-3, 'AbsTol', 1e-3, 'NMax', 20000);
    % options = bvpset('RelTol', 1e-2, 'AbsTol', 100, 'NMax', 20000);
    % xmesh = zz_wgrid;
    % solinit = bvpinit(xmesh, [0 0]);
    % p0_guess = p0;
    % solinit.y(1,:) = p0_guess;
    % solinit.y(2,2:end-1) = (p0_guess(3:end)-p0_guess(1:end-2))/dz/2; %%% Centered difference
    % solinit.y(2,1)=(p0_guess(2)-p0_guess(1))/dz;
    % solinit.y(2,end)=(p0_guess(end)-p0_guess(end-1))/dz;

    if(~hydrostatic)
        %%% non-hydrostctic
        sol1 = bvp4c(@(z,y)bvpfun(z,y,kx,zeta0), @bcfun, solinit,options);
    else
        %%% hydrostctic
        sol1 = bvp4c(@(z,y)bvpfun_hydro(z,y,zeta0), @bcfun, solinit,options);
    end
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

    U = cos(omega*t0)*Atide;
    U_wgrid = cos(omega*t0)*Atide_wgrid;


    Uzz = zeros(1,Nr+1);                     %%% on w-grid
    Uzzz = zeros(1,Nr+1);                    %%% on w-grid
    Bzz = zeros(1,Nr);                       %%% on b-grid

    if(useTanhShear)
        Uz = (U_wgrid(2:Nr+1)-U_wgrid(1:Nr))/dz; %%% on b-grid
        Uzz(2:Nr) = (Uz(2:Nr)-Uz(1:Nr-1))/dz;    %%% on w-grid
        Uzz(1)=2*Uzz(2)-Uzz(3); % Linear extrapolation
        Uzz(Nr+1)=2*Uzz(Nr)-Uzz(Nr-1);
       
        Uzzz(2:Nr)=(Uzz(3:Nr+1)-Uzz(1:Nr-1))/2/dz;    %%% on w-grid
        Uzzz(1)=2*Uzzz(2)-Uzzz(3); % Linear extrapolation
        Uzzz(Nr+1)=2*Uzzz(Nr)-Uzzz(Nr-1);
    
        Az = zeros(1,Nr+1);                            %%% on w-grid
        Az(2:Nr) = (Atide(2:Nr)-Atide(1:Nr-1))/dz;     %%% on w-grid
        Az(1)=2*Az(2)-Az(3); % Linear extrapolation
        Az(Nr+1)=2*Az(Nr)-Az(Nr-1);
        Azz = (Az(2:Nr+1)-Az(1:Nr))/dz;                %%% on b-grid
        Bzz = -Azz/omega*N^2*sind(topo)*sin(omega*t0); %%% on b-grid
    end

    dbdz(2:Nr) = (b0(2:Nr)-b0(1:Nr-1))/dz; %%% on w-grid

    %%%%%%%%%%%% B.C.-8 %%%%%%%%%%%%
    % %%% No buoyancy flux, Ocean bottom, when diffusivity is not zero
    % %%% dbdz+dB/dz+dB0/dz = 0 ==>  dbdz = -(dB/dz+dB0/dz)
    % dbdz(1) = - (omega/Shear) * (delta/Hshear) * (cosd(topo)/sind(topo)) ...
    %     +(delta/Hshear)*sin(t0); 
    % %%% No buoyancy flux, Upper boundary, when diffusivity is not zero
    % dbdz(Nr+1) = - (omega/Shear) * (delta/Hshear) * (cosd(topo)/sind(topo)) ...
    %     +(delta/Hshear)*sin(t0); 
    dbdz(1) = 0;
    dbdz(Nr+1) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    b0_wgrid(2:Nr) = 0.5*(b0(2:Nr)+b0(1:Nr-1));

    %%%%%%%%%%%% B.C.-9 %%%%%%%%%%%%
    b0_wgrid(1) = 2*b0(1)-b0_wgrid(2);
    b0_wgrid(Nr+1) = 2*b0(Nr)-b0_wgrid(Nr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    d2bdz2 = (dbdz(2:Nr+1)-dbdz(1:Nr))/dz; %%% 2nd-order centered difference

    p0_ugrid = 0.5*(p0(1:Nr)+p0(2:Nr+1));
    dpsidz = (p0(2:Nr+1)-p0(1:Nr))/dz;

    bq1(o,:) = -1i*kx*U.*b0;
    bq2(o,:) = -1i*kx*p0_ugrid*N^2*cosd(topo);
    bq3(o,:) = dpsidz*N^2*sind(topo);
    bq4(o,:) = 1i*kx*p0_ugrid.*dAdz/omega*N^2*sind(topo)*sin(omega*t0);
    bq5(o,:) = kappa*(d2bdz2-kx^2.*b0+Bzz);
 
    if(noBQ2)
        bq2(o,:) = 0;
    end
    if(noBQ3)
        bq3(o,:) = 0;
    end
    if(noBQ4)
        bq4(o,:) = 0;
    end

    dbdt(o,:) = bq1(o,:) + bq2(o,:) + bq3(o,:) + bq4(o,:) + bq5(o,:);


    % for m = 2:Nr
    %     d2zetadz2(m) = (z0(m-1)-2*z0(m)+z0(m+1))/dz^2;
    % end
    d2zetadz2(2:Nr) = (z0(1:Nr-1)-2*z0(2:Nr)+z0(3:Nr+1))/dz^2;

    %%%%%%%%%%%% B.C.-11 %%%%%%%%%%%%
    %%% Linear extrapolation
    d2zetadz2 = interp1(zz_wgrid(2:Nr),d2zetadz2(2:Nr),zz_wgrid,'linear','extrap');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    zq1(o,:) = -1i*kx*U_wgrid.*z0 + 1i*kx*p0.*Uzz;
    zq2(o,:) = 1i*kx*b0_wgrid*cosd(topo);
    zq3(o,:) = -dbdz*sind(topo);
    zq4(o,:) = nu*(d2zetadz2-kx^2.*z0-Uzzz);

    if(noZQ2)
        zq2(o,:) = 0;
    end
    if(noZQ3)
        zq3(o,:) = 0;
    end

    dzetadt(o,:) = zq1(o,:) + zq2(o,:) + zq3(o,:) + zq4(o,:);

end


%%% Code equations: see https://www.mathworks.com/help/matlab/math/solve-bvp-with-two-solutions.html
%%% Non-hydrostatic??? (why did you call it non-hydrostatic???)
%%% Re-write the relation between psi and zeta in two equations
%%% psi = psi1
%%% d(psi1)/dz = psi2
%%% d(psi2)/dz = d^2(psi)/dz^2 = kx^2*psi+zeta = kx^2*psi1+zeta
function dydz = bvpfun(z,y,kx,zeta)
    dydz = [y(2)
        kx^2*y(1)+zeta(z)];
end


%%% Hydrostatic
function dydz = bvpfun_hydro(z,y,zeta)
    dydz = [y(2)
          zeta(z)];
end


%%%%%%%%%%%% B.C.-12 %%%%%%%%%%%%
%%% Code boundary conditions for the streamfunction: phi = 0 at z=0 and z=1
function res = bcfun(ya,yb)
    res = [ya(1)
           yb(1)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    


function y = guess(x,Hmax)
   y = [sin(x/Hmax*pi)
        pi/Hmax*cos(x/Hmax*pi)];

    if(x==Hmax)
    y(1,end)=0;
    end

end