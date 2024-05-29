%%%% Integrate b and zeta with time
%%%% 1st-order and 2nd-order centered difference
%%%% Fourth-order Runge-Kutta or Third-order Adams-Bashforth
%%%% or Euler forward scheme for time advancement

clear all;
close all;

exppath = 'exps_tanh_ZeroBottom_dz2/';
constants;

% for Nexp_lambda =1:length(lambda_parm)
parfor Nexp_lambda =10:11

    lambda = lambda_parm(Nexp_lambda);
    kx = 2*pi/lambda;
    expfolder = [exppath 'lambda' num2str(lambda) '/']
    mkdir(expfolder); 

    % for Nexp_shear =1:length(Shear_parm)
    for Nexp_shear =1

        Shear = Shear_parm(Nexp_shear)
        
        expdir = [expfolder 'topo' num2str(topo) '_H' num2str(Hmax) ...
            '_N' num2str(N) '_S' num2str(Shear) ...
            '_lambda' num2str(lambda) '/'];
        outputname = [expdir 'output.mat'];
        mkdir(expdir);

        useRK4 = true;        %%% Use Tourth-order Runge-Kutta method
        useAB3 = false;       %%% Use Third-order Adams-Bashforth method

        noBQ2 = false;
        noBQ3 = false;
        noBQ4 = false;
        noZQ2 = false;
        noZQ3 = false;

        hydrostatic = false;

        if(useLinearShear)
            Atide = Shear*zz;
            Atide_wgrid = Shear*zz_wgrid;
            Umax = h_shear * Shear;
        end
        if(useTanhShear)
            %%% Zero velocity at sea floor
            Atide = h_shear*Shear *(1+ tanh( (zz  -Hmax/2) / (h_shear/2) )) /2;
            Atide_wgrid = h_shear*Shear *(1+ tanh( (zz_wgrid  -Hmax/2) / (h_shear/2) )) /2;
            Umax = h_shear * Shear;
            %%% Zero velocity at center
            % Atide = h_shear*Shear *(tanh( (zz  -Hmax/2) / (h_shear/2) )) /2;
            % Atide_wgrid = h_shear*Shear *(tanh( (zz_wgrid  -Hmax/2) / (h_shear/2) )) /2;
            % Umax = h_shear * Shear /2;
        end
        dAdz = diff(Atide_wgrid)/dz;


        %%% Background tidal velocity
        Utide =cos(tt*omega)'.*Atide;
        % Utide = repmat(cos(tt*omega)',[1 length(Atide)])...
        %     .*repmat(Atide,[length(tt) 1])/U0;
        
        
        %%%%%%%%%%%% B.C.-1 %%%%%%%%%%%%
        zeta(1,1) = 0; zeta(1,Nr+1) = 0; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for o=1:Nt-1
        
            t0 = tt(o);
            b0 = buoy(o,:);
            z0 = zeta(o,:);
            
            loop;
            psi(o,:) = p0;
        
            if(useRK4)
                RK4;
            elseif(useAB3)
                %%% Third-order Adams-Bashforth method %%%
                if (o <= 2)
                    %%% Use RK4 for the first 2 time steps
                    RK4;
                    %%% Use Euler-forward for the first 2 time steps
                    % buoy(o+1,:) = buoy(o,:) + dbdt(o,:)*dt;
                    % zeta(o+1,:) = zeta(o,:) + dzetadt(o,:)*dt;
                else
                    buoy(o+1,:) = buoy(o,:) + dt*( (23/12)*dbdt(o,:)    - (16/12)*dbdt(o-1,:)    + (5/12)*dbdt(o-2,:) );
                    zeta(o+1,:) = zeta(o,:) + dt*( (23/12)*dzetadt(o,:) - (16/12)*dzetadt(o-1,:) + (5/12)*dzetadt(o-2,:) );
                end
            else
                %%% Euler forward %%%
                buoy(o+1,:) = dbdt(o,:)*dt;
                zeta(o+1,:) = dzetadt(o,:)*dt;
            end
        
            %%%%%%%%%%%% B.C.-2 %%%%%%%%%%%%
            %%% No-stress (free-slip)
            zeta(o+1,1) = 0; zeta(o+1,Nr+1) = 0; 
        
            %%% No total stress at the ocean bottom
            % zeta(o,1) = cos(t0*omega);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            if(rem(o,round(Nt/10))==0)
                Progress = o/Nt
            end
        
        end
        
        re_psi = real(psi);   re_psi(re_psi==0)=NaN;
        re_zeta = real(zeta); re_zeta(re_zeta==0)=NaN;
        re_buoy = real(buoy); re_buoy(re_buoy==0)=NaN;
        re_dbdz = real(dbdz);
        
        re_bq1 = real(bq1);
        re_bq2 = real(bq2);
        re_bq3 = real(bq3);
        re_bq4 = real(bq4);
        re_bq5 = real(bq5);
        
        re_zq1 = real(zq1);
        re_zq2 = real(zq2);
        re_zq3 = real(zq3);
        re_zq4 = real(zq4);
        
        dbuoydz = zeros(Nt,Nr);
        for m = 2:Nr-1
            dbuoydz(:,m) = (re_buoy(:,m+1)-re_buoy(:,m-1))/dz;
        end
        
        plot_tidx = 1:10:Nt;
                
        DIV = 1;
        uuu = -real((psi(:,2:Nr+1)-psi(:,1:Nr))/dz);
        www = real(1i*kx*psi);

        zidx_b = 1:Nr;
        zidx_q = 1:Nr+1;

        bq1_int = real((sum(bq1(:,zidx_b),2)))';
        bq2_int = real((sum(bq2(:,zidx_b),2)))';
        bq3_int = real((sum(bq3(:,zidx_b),2)))';
        bq4_int = real((sum(bq4(:,zidx_b),2)))';
        bq5_int = real((sum(bq5(:,zidx_b),2)))';
        
        zq1_int = real((sum(zq1(:,zidx_q),2)))';
        zq2_int = real((sum(zq2(:,zidx_q),2)))';
        zq3_int = real((sum(zq3(:,zidx_q),2)))';
        zq4_int = real((sum(zq4(:,zidx_q),2)))';


        fit_span = round(Nt/NTtide*2)+1:Nt-1;
        
        % clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S 
        TKE = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
        TPE = 0;
        KE_PE = TKE+TPE;
        
        KE_PE_zavg = mean(KE_PE,2)';
        xxplot = tt/t1hour;
        yyplot = log(KE_PE_zavg)/2;
        [pKE,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
        % GrowthRate_KE(Nexp_lambda,Nexp_shear) = pKE(1);
        pKE(1);
        [y_fit,delta_fit] = polyval(pKE,xxplot,S);
        
        b2 = mean(re_buoy.^2,2)';
        yyplot_b2 = log(b2)/2;
        [pb2,S_b2] = polyfit(xxplot(fit_span),yyplot_b2(fit_span),1); 
        % GrowthRate(Nexp_lambda,Nexp_shear) = pb2(1);
        [y_fit_b2,delta_fit_b2] = polyval(pb2,xxplot,S_b2);

        %%% Save outputs
        s = struct('buoy',buoy,'zeta',zeta,'psi',psi, ...
            're_buoy',re_buoy,'re_psi',re_psi,'pKE',pKE,'pb2',pb2,...
            'uuu',uuu,'www',www,'NTtide',NTtide,'Nr',Nr,'Nt',Nt,'Utide',Utide,...
            'tt',tt,'t1hour',t1hour,'zz',zz,'dz',dz,'nu',nu,'kappa',kappa,'fit_span',fit_span,...
            'zidx_b',zidx_b,'bq1_int',bq1_int,'bq2_int',bq2_int,'bq3_int',bq3_int,'bq4_int',bq4_int,'bq5_int',bq5_int,...
            'zidx_q',zidx_q,'zq1_int',zq1_int,'zq2_int',zq2_int,'zq3_int',zq3_int,'zq4_int',zq4_int,...
            'Shear',Shear,'lambda',lambda,'topo',topo,'Atide',Atide,'dAdz',dAdz...
            );  
        save(sprintf(outputname),"-fromstruct",s);
        

    end
    
end


