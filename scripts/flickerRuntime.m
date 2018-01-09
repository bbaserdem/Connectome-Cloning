clear variables
file_prefix = 'fli';
save_var = {...
    'N', ...
    'D', ...
    'R', ...
    'failure', ...
    'N_values', ...
    'B_values', ...
    'D_values', ...
    'C_values', ...
    'T_values' };

N = ceil( logspace( log10(50), 3, 15) );
D = logspace( log10(0.03), log10(0.2), 8);
R = 20;

N_values = [];
B_values = [];
C_values = [];
T_values = [];

failure = zeros( length(N), length(D) );

for nN = 1:length(N)
    n = N(nN);
    for dN = 1:length(D)
        tic;
        d = D(dN);
        fprintf('Case: N=%d, D=%.2d\n', n, d);
        A = zeros(R,1);
        B = zeros(R,1);
        C = zeros(R,1);
        parfor r = 1:R
            [A(r), B(r), C(r)] = cc_flicker(n,d);
        end
        ind = A~=0;
        failure(nN,dN) = sum( ~ind ) / R;
        N_values = [N_values; n*ones(sum(ind),1) ];
        B_values = [B_values; B(ind) ];
        C_values = [C_values; A(ind) ];
        T_values = [T_values; C(ind) ];
        toc;
        fprintf( [ 'Time; ', datestr(now), '\n'] );
    end
end
D_values = B_values ./ ( N_values .^ 2 );

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