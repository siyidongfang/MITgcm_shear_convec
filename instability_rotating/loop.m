    %%% Start the loop
function [grow,vvel,buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot,pp,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert]...
        =loop(grow,j,NTtide,kappa_const,dt,Nt,dvdt,dbdt,dzetadt,omega,m0,rs,kx,shear,ss,cs,N,kappa,nu,f,tt,buoy,vvel,zeta,Diffusion,ConvectiveAdjustment,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert)


    for o=1:Nt-1

        %%% Fourth-order Runge-Kutta method %%%
        t0 = tt(o);
        b0 = buoy(o);
        z0 = zeta(o);
        v0 = vvel(o);
        [dbdt,dzetadt,dvdt]=tendency(o,dbdt,dzetadt,dvdt,omega,t0,m0,rs,kx,shear,ss,cs,N,z0,b0,v0,kappa,nu,f);
        k_1b = dbdt(o);
        k_1z = dzetadt(o);
        k_1v = dvdt(o);

        % Euler forward predictor advancing dt/2:
        b_2 = buoy(o)+0.5*dt*k_1b;
        z_2 = zeta(o)+0.5*dt*k_1z;
        v_2 = vvel(o)+0.5*dt*k_1v;
        t0 = tt(o)+dt/2;
        b0 = b_2;
        z0 = z_2;
        v0 = v_2;
        [dbdt,dzetadt,dvdt]=tendency(o,dbdt,dzetadt,dvdt,omega,t0,m0,rs,kx,shear,ss,cs,N,z0,b0,v0,kappa,nu,f);
        k_2b = dbdt(o);
        k_2z = dzetadt(o);
        k_2v = dvdt(o);

        % Euler backward corrector advancing dt/2:
        b_3 = buoy(o)+0.5*dt*k_2b;
        z_3 = zeta(o)+0.5*dt*k_2z;
        v_3 = vvel(o)+0.5*dt*k_2v;
        t0 = tt(o)+dt/2;
        b0 = b_3;
        z0 = z_3;
        v0 = v_3;
        [dbdt,dzetadt,dvdt]=tendency(o,dbdt,dzetadt,dvdt,omega,t0,m0,rs,kx,shear,ss,cs,N,z0,b0,v0,kappa,nu,f);
        k_3b = dbdt(o);
        k_3z = dzetadt(o);
        k_3v = dvdt(o);

        % Mid-point predictor advancing dt:
        b_4 = buoy(o)+dt*k_3b;
        z_4 = zeta(o)+dt*k_3z;
        v_4 = vvel(o)+dt*k_3v;
        t0 = tt(o)+dt;
        b0 = b_4;
        z0 = z_4;
        v0 = v_4;
        [dbdt,dzetadt,dvdt]=tendency(o,dbdt,dzetadt,dvdt,omega,t0,m0,rs,kx,shear,ss,cs,N,z0,b0,v0,kappa,nu,f);        
        k_4b = dbdt(o);
        k_4z = dzetadt(o);
        k_4v = dvdt(o);
    
        % Simpson rule corrector advancing dt:
        buoy(o+1) = buoy(o) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
        zeta(o+1) = zeta(o) + (1/6)*(k_1z+2*k_2z+2*k_3z+k_4z)*dt;
        vvel(o+1) = vvel(o) + (1/6)*(k_1v+2*k_2v+2*k_3v+k_4v)*dt;

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
        else
            dbdz_vert = [];
            dBdz_vert = [];
            dB0dz_vert=[];
            dbtotaldz_vert=[];
        end

    end

   
    ct = cos(omega*tt);
    st = sin(omega*tt);
    mz_t = m0-rs*st*kx;
    if(omega==0)
        mz_t = m0-shear*tt*kx;
    end
    angle_front = atand(mz_t/kx);
    a1_t = -(kx^2+mz_t.^2);  % a1_t = -(kx^2+m0^2+kx^2*rs^2*st.^2)+2*kx.*m0*rs.*st;
    psi = zeta./a1_t;
    www = 1i*kx*psi;
    uuu = -1i*mz_t.*psi;   % uuu = -1i*m0*psi+1i*kx*psi*rs.*st;

    re_buoy = real(buoy);
    re_uuu = real(uuu);
    re_www = real(www);
    % re_vvv = real(vvel);
    % pe = re_buoy.^2;
    % ke = 0.5*(re_uuu.^2+re_www.^2+re_vvv.^2);
    % kew = 0.5*(re_www.^2);
    % %%% To match Radko (2019) Eq.(19)
    if(nu*kappa~=0)
        Pr = nu/kappa;
    else
        Pr = 10;
    end

    ke_nond = (kx^2+mz_t.^2).*abs(psi/kappa_const).^2/4;    %%% To match Radko (2019) Eq.(19) 
    grav = 10;
    pe_nond = Pr*(abs(buoy)/grav).^2/4; %%% Non-dimensionalized KE and PE

    fit_span = Nt/NTtide*10+1:Nt;
    if(ConvectiveAdjustment)
        fit_span = Nt/NTtide*3+1:Nt/NTtide*10;
    end
    if(omega==0)
        fit_span = round(Nt/10):Nt;
    end
    xxplot = tt/3600;
    yyplot = log(ke_nond+pe_nond)/2;
    % yyplot = log(pe/median(pe)+ke/median(ke))/2;
    % yyplot = log(pe+ke)/2;
    % yyplot = log(pe_nond)/2;
    % yyplot = log(ke_nond)/2;
    [pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
    grow(j)=pp(1);
    % pp(1)
    % if(isnan(pp(1)))
    %     warning('NaN in growth rate!')
    % end
    % 
    % [y_fit,delta_fit] = polyval(pp,xxplot,S);
    % fig = figure(20);
    % clf;set(gcf,'Color','w')
    % plot(xxplot/24,yyplot,'LineWidth',2)
    % hold on;grid on;grid minor;
    % plot(xxplot(fit_span)/24, y_fit(fit_span));
    % hold off;
    % xlabel('Time (days)')
    % set(gca,'Fontsize',20)

    % saveas(fig,[[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '.jpeg']);


end