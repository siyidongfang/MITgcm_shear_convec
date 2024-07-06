%%%
%%% Compute velocity shear, given a Richardson number
%%%

clear


dt = 600;
NTtide = 100;
% omega = 1e-4;
omega = 2*pi/43200
Nt = NTtide/omega/dt;

topo=4;
N2 = 1e-6;

shear_convec = cosd(topo)/sind(topo)*omega;
if(topo==0)
    shear_convec = 5e-3;
end
shear_all = 0:0.25e-5:shear_convec;
ns = length(shear_all);

%%
dt_ri = dt/1000;
tt_ri = dt_ri:dt_ri:Nt*dt;
isConvec = zeros(1,ns);
Ri_min = zeros(1,ns);
for i=1:ns 
    shear = shear_all(i);
    Ri_inverse = (shear*cos(omega*tt_ri)).^2./(N2*cosd(topo) - N2*sind(topo)/omega*shear*sin(omega*tt_ri));
    % figure(50)
    % plot(tt_ri,Ri_inverse)
    % if(min(Ri_inverse)<0)
    %     isConvec(i) = 1;
    % end
    Ri_min(i) = 1/max(Ri_inverse);  
end
% Ri_min(isConvec==1)=NaN;

% plot(shear_all,1./Ri_min)


[a idx] = min(abs(Ri_min-5))
Ri5 = Ri_min(idx);
shear_Ri5 = shear_all(idx)


[a idx] = min(abs(Ri_min-1))
Ri1 = Ri_min(idx);
shear_Ri1 = shear_all(idx)


[a idx] = min(abs(Ri_min-0.25));
Ri0_25 = Ri_min(idx);
shear_Ri0_25 = shear_all(idx)

