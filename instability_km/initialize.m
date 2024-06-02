
function [Nt,tt,psi,zeta,buoy,dbdt,dzetadt,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert] = initialize(NTtide,Ptide,nt_percycle,omega,ConvectiveAdjustment,b0)

    Lt = NTtide*Ptide;
    dt = Ptide/nt_percycle;
    
    Nt = round(Lt/dt);
    tt = dt:dt:Nt*dt;
    
    if(omega==0)
        Nt = 1e3;
        dt = NTtide*Ptide/Nt;
        tt = dt:dt:Nt*dt;
    end

    psi = zeros(1,Nt);
    zeta = zeros(1,Nt);
    buoy = zeros(1,Nt);
    dbdt = zeros(1,Nt);
    dzetadt = zeros(1,Nt);
    
    if(ConvectiveAdjustment)
        dbdz_vert = zeros(1,Nt);
        dBdz_vert = zeros(1,Nt);
        dB0dz_vert = zeros(1,Nt);
        dbtotaldz_vert = zeros(1,Nt);
    else
        dbdz_vert =[];
        dBdz_vert =[];
        dB0dz_vert=[];
        dbtotaldz_vert=[];
    end
    
    %%% Initial condition
    buoy(1) = b0;
    psi(1) = 0;
    zeta(1) = 0;
    

end
