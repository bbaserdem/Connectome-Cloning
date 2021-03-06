clear variables
file_prefix = 'par';
save_var = {...
    'R', ...
    'N_values', ...
    'D_values', ...
    'C_values' };

N = ceil( logspace( 1.25, 2, 15) );
D = logspace( log10(0.03), log10(0.3), 4);
R = 20;

N_values = [];
D_values = [];
C_values = [];

for nN = 1:length(N)
    n = N(nN);
    for dN = 1:length(D)
        tic;
        d = D(dN);
        fprintf('Case: N=%d, D=%.2d\n', n, d);
        S = zeros(R,1);
        F = zeros(R,1);
        parfor r = 1:R
            [S(r), F(r)] = cc_parallel(n,d);
        end
        ind = S~=0;
        N_values = [N_values; n*ones(sum(ind),1) ];
        D_values = [D_values; F(ind) ];
        C_values = [C_values; S(ind) ];
        toc;
    end
end

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