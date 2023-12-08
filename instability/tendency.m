    


    % p0 = cumsum(cumsum(z0*dz)*dz);
    % p0 = cumsum(cumsum((kx^2*p0+z0)*dz)*dz);

    p0 = zeros(1,Nr);
    Dn = z0(2:Nr-1)'*dz^2;
    An(1,1)=C;An(1,2)=1;
    An(Nr-2,Nr-3)=1;An(Nr-2,Nr-2)=C;
    for n=2:Nr-3
        An(n,n-1)=1;
        An(n,n)=C;
        An(n,n+1)=1;
    end
    det_An = det(An);
    for n=1:Nr-2
        Bn = An;
        Bn(:,n) = Dn;
        det_Bn = det(Bn);
        p0(n+1) = det_Bn/det_An;
    end

    U = cos(t0)*Atide/U0;

    %%% Boundary condition:
    % p0(1)=0;
    % p0(Nr)=0;     
    b0(1) = b0(2); 
    b0(Nr) = b0(Nr-1); 
    % z0(1) = 3/dz^2*p0(1)-1/2*z0(2); % Woods (1954) boundary condition

    for m = 2:Nr-1
        dUdz(m) = (U(m+1)-U(m-1))/2/dz;
        dpsidz(m)   = (p0(m+1)-p0(m-1))/2/dz;
        % d2psidz2(m) = (p0(m-1)-2*p0(m)+p0(m+1))/dz^2;
        dbdz(m)   = (b0(m+1)-b0(m-1))/2/dz;
        d2bdz2(m) = (b0(m-1)-2*b0(m)+b0(m+1))/dz^2;
        d2zetadz2(m) = (z0(m-1)-2*z0(m)+z0(m+1))/dz^2;
    end

    bq1(o,:) = -1i*kx*C1*U.*b0;
    bq2(o,:) = -1i*kx*cotd(topo)*p0;
    bq3(o,:) = dpsidz;
    bq4(o,:) = +1i*kx*C1*dUdz.*tan(t0).*p0;
    bq5(o,:) = +C4*d2bdz2-C4*kx^2.*b0;
    % bq5(o,:) = zeros(1,Nr); %%% ignore diffusion

    dbdt(o,:) = bq1(o,:) + bq2(o,:) + bq3(o,:) ...
              + bq4(o,:) + bq5(o,:);

    zq1(o,:) = -1i*kx*C1*U.*z0;
    zq2(o,:) = +C2^2*(1i*kx*cotd(topo)*b0);
    zq3(o,:) = +C2^2*(-dbdz);
    zq4(o,:) = +C3*d2zetadz2-C3*kx^2.*z0;
    % zq4(o,:) = zeros(1,Nr); %%% ignore dissipation

    dzetadt(o,:) = zq1(o,:) + zq2(o,:) + zq3(o,:) + zq4(o,:);


