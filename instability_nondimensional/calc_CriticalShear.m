%%%
%%% data_processing.m
%%%

close all;clear;

topo_parm = 0:0.5:20;
Nsquare_parm = 1e1.^(-9:0.1:-4);
Shear_parm = [0.01:0.01:10 11:1:100 110:10:1000]*1e-3; 
Ptide_parm = [1 2]*43200;

NTtide = 10;
dt = 0.001;
Lt = NTtide*2*pi; % dimensionless simulation time
Nt = round(Lt/dt);
dt_real = NTtide*43200/Nt;
tt = dt:dt:Nt*dt;

NP = length(Ptide_parm);
NN = length(Nsquare_parm);
NT = length(topo_parm);
NS = length(Shear_parm);

Ri_min = NaN*zeros(NP,NN,NT,NS);
isConvec = zeros(NP,NN,NT,NS);
criticalShear = NaN*zeros(NP,NN,NT);

%%
for np = 1:NP
    np
    Ptide = Ptide_parm(np);
    omega = 2*pi/Ptide;
    for nn = 1:NN
        nn
        Nsquare = Nsquare_parm(nn);
        for nt = 1:NT
            topo = topo_parm(nt);    
            for ns=1:NS
            
                Shear = Shear_parm(ns);
            
                Ri_inverse = (Shear*cos(tt)).^2./(Nsquare*cosd(topo) - Nsquare*sind(topo)/omega*Shear*sin(tt));
                
                if(min(Ri_inverse)<0)
                    isConvec(np,nn,nt,ns) = 1;
                end
                Ri_min(np,nn,nt,ns) = 1/max(Ri_inverse);
            
            end
        end
    end
end

Ri_min(Ri_min==Inf)=NaN;


%%% Find the critical shear of Ri<=0.25
[Ri_min_critical,criticalShearidx] = min(abs(Ri_min-0.25),[],4);
for np=1:NP
    for nn=1:NN
        for nt=1:NT
            criticalShear(np,nn,nt) = Shear_parm(criticalShearidx(np,nn,nt));
        end
    end
end

criticalShear(Ri_min_critical>0.25)=NaN;
criticalShear_M2 = squeeze(criticalShear(1,:,:));
criticalShear_S1 = squeeze(criticalShear(2,:,:));

f0 = 1e-4;
Bu = NaN*zeros(NN,NT);
for nn=1:NN
    for nt=1:NT
        Bu(nn,nt) = (topo_parm(nt)/180*pi) * sqrt(Nsquare_parm(nn)) /f0; %%% Slope Burger Number
    end
end

save('CriticalShear.mat')







