

    U = cos(t0)*Atide/U0;
    p0 = cumsum(cumsum(z0*dz)*dz)-U(1)*Lz;
    p0(1)=0;
    p0(Nr)=0;     
    b0(1) = b0(2); 
    b0(Nr) = b0(Nr-1); 
    
    for m = 2:Nr-1
        dUdz(m) = (U(m+1)-U(m-1))/2/dz;
        dpsidz(m)   = (p0(m+1)-p0(m-1))/2/dz;
        d2psidz2(m) = (p0(m-1)-2*p0(m)+p0(m+1))/dz^2;
        dbdz(m)   = (b0(m+1)-b0(m-1))/2/dz;
        d2bdz2(m) = (b0(m-1)-2*b0(m)+b0(m+1))/dz^2;
        d2zetadz2(m) = (z0(m-1)-2*z0(m)+z0(m+1))/dz^2;
    end
    dbdt = dpsidz -1i*kx*cotd(topo)*p0...
         +1i*kx*Re/2*dUdz.*tan(t0).*p0 ...
         -1i*kx*Re/2*U.*b0 ...
         +(d2bdz2 - kx^2.*b0)/2/Pr; %%% dissipation
    dzetadt = 1i*kx*Re/2*dUdz.*dpsidz ...
         -1i*kx*Re/2*U.*(d2psidz2-kx^2*p0) ...
         +C^2*(1i*kx*cotd(topo)*b0-dbdz) ...
         +(d2zetadz2 - kx^2.*z0)/2; %%% dissipation


