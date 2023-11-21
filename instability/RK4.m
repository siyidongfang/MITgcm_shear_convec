    

%%% Euler forward predictor advancing dt/2:
    b_2 = buoy(o,:)+0.5*dt*k_1b;
    z_2 = zeta(o,:)+0.5*dt*k_1z;

    t0 = tt(o)+dt/2;
    b0 = b_2;
    z0 = z_2;

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
    end
    dbdt = dpsidz -1i*kx*cotd(topo)*p0...
         +1i*kx*Re/2*dUdz.*tan(t0).*p0 ...
         -1i*kx*Re/2*U.*b0 ...
         +(d2bdz2 - kx^2.*b0)/2/Pr; %%% dissipation
    dzetadt = 1i*kx*Re/2*dUdz.*dpsidz ...
         -1i*kx*Re/2*U.*(d2psidz2-kx^2*p0) ...
         +C^2*(1i*kx*cotd(topo)*b0-dbdz);

    k_2b = dbdt;
    k_2z = dzetadt;

    %%% Euler backward corrector advancing dt/2:
    b_3 = buoy(o,:)+0.5*dt*k_2b;
    z_3 = zeta(o,:)+0.5*dt*k_2z;

    t0 = tt(o)+dt/2;
    b0 = b_3;
    z0 = z_3;
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
    end
    dbdt = dpsidz -1i*kx*cotd(topo)*p0...
         +1i*kx*Re/2*dUdz.*tan(t0).*p0 ...
         -1i*kx*Re/2*U.*b0 ...
         +(d2bdz2 - kx^2.*b0)/2/Pr; %%% dissipation
    dzetadt = 1i*kx*Re/2*dUdz.*dpsidz ...
         -1i*kx*Re/2*U.*(d2psidz2-kx^2*p0) ...
         +C^2*(1i*kx*cotd(topo)*b0-dbdz);

    k_3b = dbdt;
    k_3z = dzetadt;


    %%% Mid-point predictor advancing dt:
    b_4 = buoy(o,:)+dt*k_3b;
    z_4 = zeta(o,:)+dt*k_3z;

    t0 = tt(o)+dt;
    b0 = b_4;
    z0 = z_4;
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
    end
    dbdt = dpsidz -1i*kx*cotd(topo)*p0...
         +1i*kx*Re/2*dUdz.*tan(t0).*p0 ...
         -1i*kx*Re/2*U.*b0 ...
         +(d2bdz2 - kx^2.*b0)/2/Pr; %%% dissipation
    dzetadt = 1i*kx*Re/2*dUdz.*dpsidz ...
         -1i*kx*Re/2*U.*(d2psidz2-kx^2*p0) ...
         +C^2*(1i*kx*cotd(topo)*b0-dbdz);

    k_4b = dbdt;
    k_4z = dzetadt;

    %%% Simpson rule corrector advancing dt:
    buoy(o+1,:) = buoy(o,:) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
    zeta(o+1,:) = zeta(o,:) + (1/6)*(k_1z+2*k_2z+2*k_3z+k_4z)*dt;
    