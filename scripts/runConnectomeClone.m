% Script to generate results up to a thousand neurons

N = ceil( logspace(1,3,7) )';
D = logspace(log10(0.05), log10(0.8), 6)';
repeat = 24;

Nval = zeros( length(N)*length(D), repeat );
Fval = zeros( size(Nval) );
Bval = zeros( size(Nval) );
Cval = zeros( size(Nval) );


for a = 1:length(N)
    for b = 1:length(D)
        tic;

        n = N(a);
        d = D(b);
        fprintf('\nDoing neuron: %d, density: %d\n', n, d)
        
        sto1 = zeros(1,repeat);
        sto2 = zeros(1,repeat);
        parfor e = 1:repeat
            [sto1(e), sto2(e)] = connectomeClone( n, d );
        end
        
        c = b + (a-1)*length(D);
        Nval(c,:) = n;
        Fval(c,:) = sto2 / (n^2);
        Bval(c,:) = sto2;
        Cval(c,:) = sto1;
        
        toc;
    end
end

clear a b c n d sto1 sto2 e

save results/run3.mat