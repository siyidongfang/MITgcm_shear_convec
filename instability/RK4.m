

    %%% Fourth-order Runge-Kutta method %%%
        k_1b = dbdt(o,:);
        k_1z = dzetadt(o,:);
        % Euler forward predictor advancing dt/2:
        b_2 = buoy(o,:)+0.5*dt*k_1b;
        z_2 = zeta(o,:)+0.5*dt*k_1z;
        % z_2(1) = delta/Hshear*cos(t0);
        z_2(1) = 0; z_2(Nr+1) = 0;
        t0 = tt(o)+dt/2;
        b0 = b_2;
        z0 = z_2;
        tendency;
        k_2b = dbdt(o,:);
        k_2z = dzetadt(o,:);
        % Euler backward corrector advancing dt/2:
        b_3 = buoy(o,:)+0.5*dt*k_2b;
        z_3 = zeta(o,:)+0.5*dt*k_2z;
        % z_3(1) = delta/Hshear*cos(t0);
        z_3(1) = 0; z_3(Nr+1) = 0;
        t0 = tt(o)+dt/2;
        b0 = b_3;
        z0 = z_3;
        tendency;
        k_3b = dbdt(o,:);
        k_3z = dzetadt(o,:);
        % Mid-point predictor advancing dt:
        b_4 = buoy(o,:)+dt*k_3b;
        z_4 = zeta(o,:)+dt*k_3z;
        % z_4(1) = delta/Hshear*cos(t0);
        z_4(1) = 0; z_4(Nr+1) = 0;
        t0 = tt(o)+dt;
        b0 = b_4;
        z0 = z_4;
        tendency;
        k_4b = dbdt(o,:);
        k_4z = dzetadt(o,:);
        % Simpson rule corrector advancing dt:
        buoy(o+1,:) = buoy(o,:) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
        zeta(o+1,:) = zeta(o,:) + (1/6)*(k_1z+2*k_2z+2*k_3z+k_4z)*dt;

