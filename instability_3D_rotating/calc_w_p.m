function [w0 p0]=calc_w_p(o,t0,m0,k0,l0,u0,v0,b0,ss,cs,shear,omega,f0,nu)

    st = sin(omega*t0);
    ct = cos(omega*t0);
    mz = m0 - shear/omega*st*k0 -shear*f0*cs/omega^2*(ct-1)*l0;
    D = -(k0^2+l0^2+mz^2);

    w0 = -(k0*u0+l0*v0)/mz;
    dwdt = (w0-w(o-1))/dt;
    p0 = (b0*cs+nu*D*w0-dwdt-f0*v0*ss)/1i/mz;

end