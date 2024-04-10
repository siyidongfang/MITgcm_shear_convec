    %%% Start the loop

    for o=1:Nt-1
        %%% Fourth-order Runge-Kutta method %%%
        t0 = tt(o);
        b0 = buoy(o);
        z0 = zeta(o);
        tendency;
        k_1b = dbdt(o);
        k_1z = dzetadt(o);
        % Euler forward predictor advancing dt/2:
        b_2 = buoy(o)+0.5*dt*k_1b;
        z_2 = zeta(o)+0.5*dt*k_1z;
        t0 = tt(o)+dt/2;
        b0 = b_2;
        z0 = z_2;
        tendency;
        k_2b = dbdt(o);
        k_2z = dzetadt(o);
        % Euler backward corrector advancing dt/2:
        b_3 = buoy(o)+0.5*dt*k_2b;
        z_3 = zeta(o)+0.5*dt*k_2z;
        t0 = tt(o)+dt/2;
        b0 = b_3;
        z0 = z_3;
        tendency;
        k_3b = dbdt(o);
        k_3z = dzetadt(o);
        % Mid-point predictor advancing dt:
        b_4 = buoy(o)+dt*k_3b;
        z_4 = zeta(o)+dt*k_3z;
        t0 = tt(o)+dt;
        b0 = b_4;
        z0 = z_4;
        tendency;
        k_4b = dbdt(o);
        k_4z = dzetadt(o);
    
        % Simpson rule corrector advancing dt:
        buoy(o+1) = buoy(o) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
        zeta(o+1) = zeta(o) + (1/6)*(k_1z+2*k_2z+2*k_3z+k_4z)*dt;

        if(ConvectiveAdjustment)
            dbdz_vert(o) = 1/cs * real(1i*mz*buoy(o));
            dBdz_vert(o) = 1/cs * (-rs*N^2*ss*st); %%% Vertical direction
            dB0dz_vert(o) =1/cs * (N^2*cs);
            dbtotaldz_vert(o) = dB0dz_vert(o)+dBdz_vert(o)+dbdz_vert(o);
            
            if(dbtotaldz_vert(o)<=0)
                nu=0.01;
                kappa=0.01;
                % imagb = 1/mz*(-rs*N^2*ss*st + N^2*cs);
                % buoy(o+1)=real(buoy(o))+1i*imagb;
            else
                if(Diffusion)
                    kappa = kappa_const;
                    nu = nu_const;
                else 
                    kappa = 0;
                    nu = 0;
                end
            end
        end

    end

    
    ct = cos(omega*tt);
    st = sin(omega*tt);
    mz_t = m0-rs*st*kx;
    a1_t = -(kx^2+mz_t.^2);  % a1_t = -(kx^2+m0^2+kx^2*rs^2*st.^2)+2*kx.*m0*rs.*st;
    psi = zeta./a1_t;
    www = 1i*kx*psi;
    uuu = -1i*mz_t.*psi;   % uuu = -1i*m0*psi+1i*kx*psi*rs.*st;

    re_buoy = real(buoy);
    re_uuu = real(uuu);
    re_www = real(www);
    pe = re_buoy.^2;
    ke = 0.5*(re_uuu.^2+re_www.^2);
    kew = 0.5*(re_www.^2);
    % %%% To match Radko (2019) Eq.(19)
    if(nu*kappa~=0)
        Pr = nu/kappa;
    else
        Pr = 1;
    end
    pe = Pr*pe/4; %%% To match Radko (2019) Eq.(19)
    ke = ke/2; 
    kew = kew/2;
    
    % fit_span = Nt/NTtide*3+1:Nt/NTtide*10;
    fit_span = Nt/NTtide*3+1:Nt;
    if(omega==0)
        fit_span = 1:Nt;
    end
    xxplot = tt/3600;
    yyplot = log(pe/median(pe)+ke/median(ke))/2;
    % yyplot = log(pe+ke)/2;
    [pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
    grow(i) = pp(1);
    if(isnan(grow(i)))
        warning('NaN in growth rate!')
    end
    % [y_fit,delta_fit] = polyval(pp,xxplot,S);
    % figure(20)
    % clf;
    % plot(xxplot,yyplot)
    % hold on;grid on;grid minor;
    % plot(xxplot(fit_span), y_fit(fit_span));
    % hold off;

