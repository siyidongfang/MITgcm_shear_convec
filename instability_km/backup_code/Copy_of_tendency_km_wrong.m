    

    %%% db/dt
    bq1(k,m,d,o) = -1i*kx(k)*C1*z(d)*cos(t0).*b0;
    bq2(k,m,d,o) = -1i*kx(k)*cotd(topo).*p0;
    bq3(k,m,d,o) = +1i*mz(m)*p0;
    bq4(k,m,d,o) = +1i*kx(k)*C1.*sin(t0).*p0;
    bq5(k,m,d,o) = -C4*(kx(k)^2+mz(m)^2).*b0;
    
    dbdt(k,m,d,o) = bq1(k,m,d,o) + bq2(k,m,d,o) + ...
        bq3(k,m,d,o) + bq4(k,m,d,o) + bq5(k,m,d,o);
    
    %%% d(psi)/dt
    pq1(k,m,d,o) = -1i*kx(k)*C1*z(d)*cos(t0).*p0;
    pq2(k,m,d,o) = -C2*1i*kx(k)*cotd(topo)/(kx(k)^2+mz(m)^2)*b0;
    pq3(k,m,d,o) = +C2*1i*mz(m)/(kx(k)^2+mz(m)^2)*b0;
    pq4(k,m,d,o) = -C3*(kx(k)^2+mz(m)^2)*p0;

    dpdt(k,m,d,o) = pq1(k,m,d,o) + pq2(k,m,d,o) + pq3(k,m,d,o) + pq4(k,m,d,o);