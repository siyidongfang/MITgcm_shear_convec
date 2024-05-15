    % % % % for m = 2:Nr-1
    % % % %     d2bdz2(m) = (b0(m-1)-2*b0(m)+b0(m+1))/dz^2;
    % % % % end
    % % % % 
    % % % % %%%%%%%%%%%% B.C.-10 %%%%%%%%%%%%
    % % % % dbdz1 = 0.5*(dbdz(1)+dbdz(2));
    % % % % dbdzNr = 0.5*(dbdz(Nr)+dbdz(Nr+1));
    % % % % 
    % % % % bw2 = b0_wgrid(2); 
    % % % % b1 = b0(1);
    % % % % zw2 = zz_wgrid(2);
    % % % % z1 = zz(1);
    % % % % 
    % % % % bwNr = b0_wgrid(Nr); 
    % % % % bNr = b0(Nr);
    % % % % zwNr = zz_wgrid(Nr);
    % % % % zNr = zz(Nr);
    % % % % 
    % % % % % Taylor expansion: bw2 ~= b1 + dbdz1*(zw2-z1) + 0.5*d2bdz2(1)*(zw2-z1)^2
    % % % % d2bdz2(1) = ( bw2 - b1 - dbdz1*(zw2-z1) ) ./ ( 0.5*(zw2-z1)^2 );
    % % % % % Taylor expansion: bwNr ~= bNr + dbdzNr*(zwNr-zNr) + 0.5*d2bdz2(Nr)*(zwNr-zNr)^2
    % % % % d2bdz2(Nr) = ( bwNr - bNr - dbdzNr*(zwNr-zNr) ) ./ ( 0.5*(zwNr-zNr)^2 );
    % % % % 
    % % % % %%% Linear extrapolation
    % % % % % d2bdz2([1 Nr]) = interp1(zz(2:Nr-1),d2bdz2(2:Nr-1),zz([1 Nr]),'linear','extrap');
    % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
