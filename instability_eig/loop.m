    %%% Start the loop
    function [buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t] ...
              =loop(dt,Nt,dbdt,dzetadt,omega,m0,rs,kx,shear,ss,cs,N,kappa,nu,tt,buoy,zeta)

    for o=1:Nt-1
        %%% Fourth-order Runge-Kutta method %%%
        t0 = tt(o);
        b0 = buoy(o);
        z0 = zeta(o);
        [dbdt,dzetadt]=tendency(o,dbdt,dzetadt,omega,t0,m0,rs,kx,shear,ss,cs,N,z0,b0,kappa,nu);
        k_1b = dbdt(o);
        k_1z = dzetadt(o);
        % Euler forward predictor advancing dt/2:
        b_2 = buoy(o)+0.5*dt*k_1b;
        z_2 = zeta(o)+0.5*dt*k_1z;
        t0 = tt(o)+dt/2;
        b0 = b_2;
        z0 = z_2;
        [dbdt,dzetadt]=tendency(o,dbdt,dzetadt,omega,t0,m0,rs,kx,shear,ss,cs,N,z0,b0,kappa,nu);
        k_2b = dbdt(o);
        k_2z = dzetadt(o);
        % Euler backward corrector advancing dt/2:
        b_3 = buoy(o)+0.5*dt*k_2b;
        z_3 = zeta(o)+0.5*dt*k_2z;
        t0 = tt(o)+dt/2;
        b0 = b_3;
        z0 = z_3;
        [dbdt,dzetadt]=tendency(o,dbdt,dzetadt,omega,t0,m0,rs,kx,shear,ss,cs,N,z0,b0,kappa,nu);
        k_3b = dbdt(o);
        k_3z = dzetadt(o);
        % Mid-point predictor advancing dt:
        b_4 = buoy(o)+dt*k_3b;
        z_4 = zeta(o)+dt*k_3z;
        t0 = tt(o)+dt;
        b0 = b_4;
        z0 = z_4;
        [dbdt,dzetadt]=tendency(o,dbdt,dzetadt,omega,t0,m0,rs,kx,shear,ss,cs,N,z0,b0,kappa,nu);        
        k_4b = dbdt(o);
        k_4z = dzetadt(o);
    
        % Simpson rule corrector advancing dt:
        buoy(o+1) = buoy(o) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
        zeta(o+1) = zeta(o) + (1/6)*(k_1z+2*k_2z+2*k_3z+k_4z)*dt;

    end

    ct = cos(omega*tt);
    st = sin(omega*tt);
    mz_t = m0-rs*st*kx;
    if(omega==0)
        mz_t = m0-shear*tt*kx;
    end
    angle_front = atand(mz_t/kx);
    a1_t = -(kx^2+mz_t.^2);  % a1_t = -(kx^2+m0^2+kx^2*rs^2*st.^2)+2*kx.*m0*rs.*st;
    psi = zeta./a1_t;
    www = 1i*kx*psi;
    uuu = -1i*mz_t.*psi;   % uuu = -1i*m0*psi+1i*kx*psi*rs.*st;

    re_buoy = real(buoy);
    re_uuu = real(uuu);
    re_www = real(www);


end