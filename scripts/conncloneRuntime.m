% Script to generate results up to a thousand neurons
clear variables
file_prefix = 'run';
save_var = {...
    'N', ...
    'D', ...
    'repeat', ...
    'Nval', ...
    'Dval', ...
    'Cval', ...
    'Bval' };

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