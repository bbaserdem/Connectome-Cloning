% Script to generate results up to a thousand neurons with scale-free
% networks
clear variables
file_prefix = 'sfn';
save_var = {...
    'N_values', ...
    'D_values', ...
    'C_values', ...
    'regRes', ...
    'regResConf', ...
    'regResStats' };

N = ceil( logspace(log10(20),3,7) )';
K = 1:10;
repeat = 24;

N_values = [];
D_values = [];
C_values = [];

for nN = 1:length(N)
    n = N(nN);
    for kN = 1:length(K)
        k = K(kN);
        fprintf('\nDoing neuron: %d, degree: %d\n', n, k);
        tic;
        sto = zeros(repeat,1);
        parfor e = 1:repeat
            sto(e) = connectomeClone( n, k );
        end
        ind = sto~=0;
        N_values = [N_values; n * ones(sum(ind),1) ];
        D_values = [D_values; (k/n) * ones(sum(ind),1) ];
        C_values = [C_values; sto(ind)];
        toc;
    end
end

[regRes, regResConf, ~, ~, regResStats] = regress( ...
    log(C_values), ...
    [ ones(size(C_values)), log(D_values), log(N_values) ] );
regRes(1) = exp(regRes(1));
regResConf(1,:) = exp(regResConf(1,:));

file = [ 'results/', file_prefix, ...
    sprintf( '%02d', 1 + length(dir(['results/',file_prefix,'*'])) ), ...
    '.mat' ];
for a = 1:length(save_var)
    if a == 1
        save( file, save_var{a} );
    else
        save( file, save_var{a}, '-append' );
    end
end
system( [ 'git add ', file ] );