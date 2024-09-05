    %%% Start the loop
function [grow,bb,uu,psi,www,uuu,re_bb,re_uu,re_ww,ct,st,mz_t,angle_front,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot,pp,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert]...
        =loop(grow,j,NTtide,kappa_const,dt,Nt,dbdt,dzetadt,omega,m0,rs,k0,shear,ss,cs,N,kappa,nu,tt,bb,uu,Diffusion,ConvectiveAdjustment,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert)


    for o=1:Nt-1
        
        %%% Fourth-order Runge-Kutta method %%%
        t0 = tt(o);
        b0 = bb(o);
        u0 = uu(o);
        v0 = vv(o);
        [dbdt,dzetadt]=tendency(o,dbdt,dzetadt,omega,t0,m0,rs,k0,shear,ss,cs,N,u0,b0,kappa,nu);
        k_1b = dbdt(o);
        k_1z = dzetadt(o);

        % Euler forward predictor advancing dt/2:
        b_2 = bb(o)+0.5*dt*k_1b;
        z_2 = uu(o)+0.5*dt*k_1z;
        t0 = tt(o)+dt/2;
        b0 = b_2;
        u0 = z_2;
        [dbdt,dzetadt]=tendency(o,dbdt,dzetadt,omega,t0,m0,rs,k0,shear,ss,cs,N,u0,b0,kappa,nu);
        k_2b = dbdt(o);
        k_2z = dzetadt(o);

        % Euler backward corrector advancing dt/2:
        b_3 = bb(o)+0.5*dt*k_2b;
        z_3 = uu(o)+0.5*dt*k_2z;
        t0 = tt(o)+dt/2;
        b0 = b_3;
        u0 = z_3;
        [dbdt,dzetadt]=tendency(o,dbdt,dzetadt,omega,t0,m0,rs,k0,shear,ss,cs,N,u0,b0,kappa,nu);
        k_3b = dbdt(o);
        k_3z = dzetadt(o);

        % Mid-point predictor advancing dt:
        b_4 = bb(o)+dt*k_3b;
        z_4 = uu(o)+dt*k_3z;
        t0 = tt(o)+dt;
        b0 = b_4;
        u0 = z_4;
        [dbdt,dzetadt]=tendency(o,dbdt,dzetadt,omega,t0,m0,rs,k0,shear,ss,cs,N,u0,b0,kappa,nu);        k_4b = dbdt(o);
        k_4z = dzetadt(o);
    
        % Simpson rule corrector advancing dt:
        bb(o+1) = bb(o) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
        uu(o+1) = uu(o) + (1/6)*(k_1z+2*k_2z+2*k_3z+k_4z)*dt;

    end

   

    re_bb = real(bb);
    re_uu = real(uu);
    re_ww = real(ww);
    re_vv = real(vv);
    re_pp = real(pp);

    % pe = re_buoy.^2;
    % ke = 0.5*(re_uuu.^2+re_www.^2);
    % %%% To match Radko (2019) Eq.(19)

    if(nu*kappa~=0)
        Pr = nu/kappa;
    else
        Pr = 1;
    end

    ke_nond = (k0^2+mz_t.^2).*abs(psi/kappa_const).^2/4;    %%% To match Radko (2019) Eq.(19) 
    grav = 10;
    pe_nond = Pr*(abs(bb)/grav).^2/4; %%% Non-dimensionalized KE and PE

    fit_span = Nt/NTtide*1+1:Nt;
    xxplot = tt/3600;
    yyplot = log(ke_nond+pe_nond)/2;
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