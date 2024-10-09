    clear;
    constants_linear
    xmesh = zz_wgrid;
    % solinit = bvpinit(xmesh, [0 0]);
    solinit = bvpinit(xmesh,@(x)guess(x,Hmax));


function y = guess(x,Hmax)
   y = [sin(x/Hmax*pi)
        pi/Hmax*cos(x/Hmax*pi)];

    if(x==Hmax)
    y(1,end)=0;
    end

end