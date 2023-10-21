%%%
%%% createRunName.m
%%%
%%% Standardizes generation of simulation names based on model parameters
%%%

function run_name = createRunName (Atide,randtopog_height,randtopog_length,Nr,Nx,run_type)

% run_name=['RH',num2str(randtopog_height),'RL',num2str(randtopog_length),...
%     '_Atide',num2str(Atide),'_Nr',num2str(Nr),'Ny',num2str(Ny),'_',run_type];
run_name=['RH',num2str(randtopog_height),'RL',num2str(randtopog_length),...
    '_Atide',num2str(Atide),'_Nr',num2str(Nr),'Nx',num2str(Nx)];
end