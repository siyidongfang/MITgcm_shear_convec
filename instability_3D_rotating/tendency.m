
function [ww,pp,dwdt,dudt,dvdt,dbdt]=tendency(o,RKstep,ww,pp,dwdt,dudt,dvdt,dbdt,t0,m0,k0,l0,u0,v0,b0,ss,cs,shear,omega,N,f0,kappa,nu)

    st = sin(omega*t0);
    ct = cos(omega*t0);
    mz = m0 - shear/omega*st*k0 -shear*f0*cs/omega^2*(ct-1)*l0;
    D = -(k0^2+l0^2+mz^2);
    Uz = shear*ct;
    Vz = - shear*f0*cs/omega*st;
    Bz = -shear/omega*N^2*ss*st;

    w0 = -(k0*u0+l0*v0)/mz;

    if(RKstep==1)
        ww(o) = w0;
        dwdt(o) = (w0-ww(o-1))/dt;
    elseif(RKstep==2)
        dwdt(o) = (w0-ww(o))/(dt/2);
    elseif(RKstep==3)
        dwdt(o) = (w0-ww(o))/(dt/2);
    elseif(RKstep==4)
        dwdt(o) = (w0-ww(o))/dt;
    end

    p0 = (b0*cs+nu*D*w0-dwdt(o)-f0*v0*ss)/1i/mz;
    if(RKstep==1)
        pp(o) = p0;
    end

    dudt(o) = - w0*Uz + f0*v0*cs - 1i*k0*p0 + b0*ss + nu*D*u0;
    dvdt(o) = - w0*Vz - f0*u0*cs - f0*w0*ss - 1i*l0*p0 + nu*D*v0;
    dbdt(o) = - u0*N^2*ss - w0*N^2*cs - w0*Bz + kappa*D*b0;

end




