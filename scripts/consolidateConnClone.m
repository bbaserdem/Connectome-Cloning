% This script consolidates data from different trials
% Consolidates Nval, Fval and Cval.

clear variables
load results/run_all.mat
load results/run3.mat Nval Fval Cval

% Add results to database
mask = Cval(:)~=0;
N_values = [ N_values; ( Nval( mask ) ) ];
D_values = [ D_values; ( Fval( mask ) ) ];
C_values = [ C_values; ( Cval( mask ) ) ];

clear Nval Fval Cval mask

% Run new regression, of the form C = X1 * D^X2 * N^X3 

[regRes, regResConf, ~, ~, regResStats] = regress( ...
    log(C_values), ...
    [ ones(size(C_values)), log(D_values), log(N_values) ] );
regRes(1) = exp(regRes(1));
regResConf(1,:) = exp(regResConf(1,:));

% Save new results
save results/run_all.mat