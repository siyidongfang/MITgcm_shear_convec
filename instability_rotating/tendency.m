
function [dbdt,dzetadt,dvdt]=tendency(o,dbdt,dzetadt,dvdt,omega,t0,m0,rs,kx,shear,ss,cs,N,z0,b0,v0,kappa,nu,f)

    st = sin(omega*t0);
    mz = m0-rs*st*kx;
    
    if(omega==0)
        mz = m0-shear*t0*kx;
    end

    a0 = (1i*m0*ss-1i*kx*cs)*N^2; 
    a1 = -(kx^2+mz^2); 
    a2 = 1i*kx*cs - 1i*mz*ss;
    a3 = -1i*kx*f*ss - 1i*mz*f*cs;
    a4 = -1i*kx*f*ss + 1i*mz*f*cs;
    
    p0 = z0/a1;
    
    dbdt(o) = a0*p0 + kappa*a1*b0;
    dzetadt(o) = a2*b0 + a3*v0 + nu*a1*z0;
    dvdt(o) = a4*p0 + nu*a1*v0;

end




