
function [dbdt,dzetadt]=tendency(o,dbdt,dzetadt,omega,t0,m0,rs,kx,shear,ss,cs,N,z0,b0,kappa,nu)

    st = sin(omega*t0);
    mz = m0-rs*st*kx;
    
    if(omega==0)
        mz = m0-shear*t0*kx;
    end

    a0 = (1i*m0*ss-1i*kx*cs)*N^2; 
    a1 = -(kx^2+mz^2);        % a1 = -(kx^2+m0^2+kx^2*rs^2*st^2)+2*kx*m0*rs*st;
    a2 = 1i*kx*cs - 1i*mz*ss; % a3 = 1i*kx*(cs+rs*ss*st) - 1i*m0*ss;
    
    p0 = z0/a1;
    
    dbdt(o) = a0*p0 + kappa*a1*b0;
    dzetadt(o) = a2*b0 + nu*a1*z0;


end




