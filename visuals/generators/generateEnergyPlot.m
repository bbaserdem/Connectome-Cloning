% USe this script to generate the energy plots.
% Will need paralelization for 1000x1000

N = ceil( logspace(1,1,1) );
D = 0.1;
T = 1000;
R = 10;

cEstimate = @(n,d) round( exp(2.5199) * ( d ^ 1.6309 ) * ( n ^ 3.3652 ) );
sEstimate = @(n,d) ceil( cEstimate(n,d) / T );
Estimate = @(n,d) length( 1:sEstimate(n,d):(2*cEstimate(n,d)) );
rec = cell(length(D), length(N), 2);

for d_n = 1:length(D)
    d = D(d_n)
    for n_n = 1:length(N)
        n = N(n_n);
        tic;
        fprintf('Doing N=%d, D=%.2d\n', n, d);
        C = Estimate(n,d);
        X = zeros(C,R);
        Y = zeros(C,R);
        parfor r = 1:R
            [ X(:,r), Y(:,r) ] = cc_energyPlot(n, d, T);
        end
        rec{d_n, n_n, 1} = X;
        rec{d_n, n_n, 2} = Y;
        toc;
    end
end

colors = hsv2rgb( [ (1:length(N))'/length(N), 0.8*ones(length(N),2) ] );
labels = ['N=', num2str(N')];

% Plot now
figure(1);
hold on;
for d_n = 1:length(D)
    for n_n = 1:length(D)
        for r = 1:R
            X = rec{d_n, n_n, 1}(rec{d_n, n_n, 1}(:,r)~=0, r);
            Y = rec{d_n, n_n, 2}(rec{d_n, n_n, 1}(:,r)~=0, r) + 1;
            plot( X, Y, 'Color', colors(n_n,:), 'LineWidth', 3 );
        end
    end
end
hold off
title('Energy vs time step');
set( gca, 'xscale', 'log', 'yscale', 'log');
legend(labels);
xlabel('Time (steps)');
ylabel('Energy');
set(gca,'FontSize',25);