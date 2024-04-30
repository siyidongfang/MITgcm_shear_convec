    


    % p0 = cumsum(cumsum(z0*dz)*dz);
    % p0 = cumsum(cumsum((kx^2*p0+z0)*dz)*dz);

    % z0(1) = 3/dz^2*p0(2)-1/2*z0(2); % Woods (1954) boundary condition for extrapolation
    % z0(1) = p0(2)/dz^2;
    % z0(Nr+1) = p0(Nr)/dz^2; 

    Dn = z0(2:Nr-1)'*dz^2;
    An(1,1)=C;An(1,2)=1;
    An(Nr-2,Nr-3)=1;An(Nr-2,Nr-2)=C;
    for n=2:Nr-3
        An(n,n-1)=1;
        An(n,n)=C;
        An(n,n+1)=1;
    end
    det_An = det(An);

    if(useParallel)
        parfor n=1:Nr-2
            Bn = An;
            Bn(:,n) = Dn;
            det_Bn = det(Bn);
            p0(n+1) = det_Bn/det_An;
        end
    else
        for n=1:Nr-2
            Bn = An;
            Bn(:,n) = Dn;
            det_Bn = det(Bn);
            p0(n+1) = det_Bn/det_An;
        end
    end

    U = cos(t0)*Atide/U0;
    U_wgrid = cos(t0)*Atide_wgrid/U0;
    if(U0==0)
        U=zeros(1,Nr);
        U_wgrid = zeros(1,Nr+1);
    end

    %%% Boundary condition:
    p0(1)=0;
    p0(Nr+1)=0;     
    b0(1) = b0(2); 
    b0(Nr) = b0(Nr-1);

    for m = 2:Nr
        % dpsidz(m)   = (p0(m+1)-p0(m-1))/2/dz;
        % d2psidz2(m) = (p0(m-1)-2*p0(m)+p0(m+1))/dz^2;
        d2zetadz2(m) = (z0(m-1)-2*z0(m)+z0(m+1))/dz^2;
    end

    for m = 2:Nr-1
        dUdz(m) = (U(m+1)-U(m-1))/2/dz;
        % dbdz(m)   = (b0(m+1)-b0(m-1))/2/dz;
        d2bdz2(m) = (b0(m-1)-2*b0(m)+b0(m+1))/dz^2;
    end

    dUdz(1) = dUdz(2); dUdz(Nr) = dUdz(Nr-1);
    dbdz(1) = 0; dbdz(Nr) = 0;
    d2bdz2(1) = 0; d2bdz2(Nr) = 0;

    % dpsidz(1) = (p0(2)-p0(1))/dz; %%% 1st order difference for boundaries
    % dpsidz(Nr+1) = (p0(Nr+1)-p0(Nr))/dz;
    d2zetadz2(1) = d2zetadz2(2);d2zetadz2(Nr+1) = d2zetadz2(Nr); %????? Is this correct?

    p0_ugrid = 0.5*(p0(1:Nr)+p0(2:Nr+1));
    dpsidz = (p0(2:Nr+1)-p0(1:Nr))/dz;

    bq1(o,:) = -1i*kx*C1*U.*b0;
    bq2(o,:) = -1i*kx*cotd(topo).*p0_ugrid;
    bq3(o,:) = dpsidz;
    bq4(o,:) = +1i*kx*C1*dUdz.*tan(t0).*p0_ugrid;
    % bq5(o,:) = +C4*d2bdz2-C4*kx^2.*b0;
    if(NOdiffusion)
        bq5(o,:) = zeros(1,Nr); %%% ignore diffusion
    end

    dbdt(o,:) = bq1(o,:) + bq2(o,:) + bq3(o,:) ...
              + bq4(o,:) + bq5(o,:);

    dbdz = zeros(1,Nr+1);
    dbdz(2:Nr) = (b0(2:Nr)-b0(1:Nr-1))/dz;
    b0_wgrid = zeros(1,Nr+1);
    b0_wgrid(2:Nr) = 0.5*(b0(2:Nr)+b0(1:Nr-1));

    zq1(o,:) = -1i*kx*C1*U_wgrid.*z0;
    zq2(o,:) = +C2^2*(1i*kx*cotd(topo)*b0_wgrid);
    zq3(o,:) = +C2^2*(-dbdz);
    % zq4(o,:) = +C3*d2zetadz2-C3*kx^2.*z0;
    if(NOdiffusion)
        zq4(o,:) = zeros(1,Nr+1); %%% ignore dissipation
    end

    dzetadt(o,:) = zq1(o,:) + zq2(o,:) + zq3(o,:) + zq4(o,:);


