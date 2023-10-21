%%%
%%% load_spacing.m
%%%


    %%% Grid spacing matrices
    DX_xy = repmat(delX',[1 Ny]);
    DY_xy = repmat(delY,[Nx 1]);
    DY_yz = repmat(delY',[1 Nr]);
    DZ_yz = repmat(delR,[Ny 1]);
    DX = repmat(reshape(delX,[Nx 1 1]),[1 Ny Nr]);
    DY = repmat(reshape(delY,[1 Ny 1]),[Nx 1 Nr]);
    DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);

    dy = delY(1);
    dx = delX(1);

    [YY,XX] = meshgrid(yy,xx);



