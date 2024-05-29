%%% 
%%% decompose.m
%%%
%%% Calculate the percentage of integrated forcing, 
%%% following Kaiser & Pratt 2022.
%%% Run this script after numerical.m

%%% To estimate the relative contribution of each term to the instability:
%%% For different wavenumber kx
%%% Vertically integrate the right-hand-side terms of dbdt, dzetadt
%%% Then cumulatively integrate those terms with time

