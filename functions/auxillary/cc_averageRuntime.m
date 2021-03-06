function [ T ] = cc_averageRuntime( N, D )
% CC_AVERAGERUNTIME Automatically generates projected runtime
% Loads results from results/run_all, a

load( [getenv('HOME'), '/Dropbox/lab_projects/bbaserde/', ...
    'ConnectomeCloning/results/cc_all.mat'], 'regRes' );
T = round( regRes(1) * (D^regRes(2)) * (N^regRes(3)) );

end

