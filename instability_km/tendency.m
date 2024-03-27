

ct = cos(omega*t0);
st = sin(omega*t0);

a1 = -(kx^2+mz^2+2*kx*mz*rs*st+kx^2*rs^2*st.^2);
a2 = 1i*kx*(2*rs*ss*st-cs)*N^2 + 1i*mz*N^2*ss;
a3 = 1i*kx*(cs-rs*ss*st) - 1i*mz*ss;

p0 = z0/a1;

dbdt(o) = a2*p0;

dzetadt(o) = a3*b0;

