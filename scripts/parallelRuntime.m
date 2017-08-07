clear variables

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

save results/par2.mat N_values D_values C_values