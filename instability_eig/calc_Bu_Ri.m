

%%% Calculate the Slope Burger number and Richardson number

load('products/grow_K1_new.mat')
constants_sens_K1

% load('products/grow_M2_new.mat')
% constants_sens_M2

% load('products/grow_ref.mat')
% constants_ref;

f0 = 1e-4;

Bu = NaN.*zeros(Ntopo,Nn,Ns);
Ri = NaN.*zeros(Ntopo,Nn,Ns);
strat = NaN.*zeros(Ntopo,Nn,Ns);
isConvec = NaN.*zeros(Ntopo,Nn,Ns);

velocityshear = NaN.*zeros(Ntopo,Nn,Ns);

dt_ri = 10;
Nt_ri = 100*Ptide/dt_ri;
tt_ri = dt_ri:dt_ri:Nt_ri*dt_ri;


for ntopo = 1:Ntopo

    ntopo
    topo = topo_all(ntopo);
    alpha = topo/180*pi; 

    for nn = 1:Nn
        N = N_all(nn);
        % Slope Burger number
        Bu(ntopo,nn,:) = N*alpha/f0;
        strat(ntopo,nn,:) = N;

        clear shear_all;
        shear_all = [0:3*N/(Ns-1):3*N]; 
        Ns = length(shear_all);

       parfor ns =1:Ns
            shear = shear_all(ns);
            % Richardson number
            Ri_inverse = (cosd(topo)*shear*cos(omega*tt_ri)).^2./(N^2 - N^2*sind(topo)*cosd(topo)/omega*shear*sin(omega*tt_ri));
            % figure(50);plot(tt_ri,Ri_inverse)
            Ri(ntopo,nn,ns) = 1/max(Ri_inverse);  
            if(min(Ri_inverse)<0)
                isConvec(ntopo,nn,ns) = 1;
            end

            velocityshear(ntopo,nn,ns) = shear;

        end

    end

end


save('products/grow_K1_calc_new.mat')


% % for ntopo = 1:Ntopo
% for ntopo = 2
%     grow_n = squeeze(grow_eig(ntopo,:,:));
%     figure(1)
%     pcolor(N_all,shear_all,grow_n')
%     shading flat;
%     colorbar;
% end