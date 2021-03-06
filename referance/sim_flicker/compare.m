Ndensity_con       = cell( length(N), length(D) );
Nconvergence       = cell( length(N), length(D) );
[Ndensity_con{:}]  = deal( zeros(1,     repeat) );
[Nconvergence{:}]  = deal( zeros(1,     repeat) );

fprintf('Doing regular network\n');
for n = 1:length(N)
    nt = N(n);
    fprintf('NEURONS NO:\t%d\n',nt);
    for d = 1:length(D)
        dt = D(d);
        fprintf('Density:%.1d\n',dt);
        sto1 = zeros(1,repeat);
        sto2 = zeros(1,repeat);
        tic;
        for r = 1:repeat
            [Nconvergence{n,d}(r), Ndensity_con{n,d}(r)] = ...
                sim_erdos( nt, dt );
        end
        toc;
    end
end

Ncolors = hsv2rgb( [...
    (1:length(D))'/length(D), 0.75*ones(length(D),2) ] );
labels = cell(1, 2*length(D) );

figure(6);
hold on;
for d = 1:length(D)
    labels{d} = ['(jitter) ', num2str(D(d))];
    x = [];
    y = [];
    for n = 1:length( N )
        x = [ x, (n^4) * density_con{n,d}(convergence{n,d}~=0)];
        y = [ y,         convergence{n,d}(convergence{n,d}~=0)];
    end
    scatter(...
        x,...
        y,...
        10,...
        colors(d,:), 'filled' ...
        );
end
for d = 1:length(D)
    labels{length(D)+d} = num2str(D(d));
    x = [];
    y = [];
    for n = 1:length( N )
        x = [ x, (n^4) * Ndensity_con{n,d}(Nconvergence{n,d}~=0)];
        y = [ y,         Nconvergence{n,d}(Nconvergence{n,d}~=0)];
    end
    scatter(...
        x,...
        y,...
        10,...
        Ncolors(d,:), 'filled' ...
        );
end
set( gca, 'xscale', 'log', 'yscale', 'log');
title('Convergence (Flickering synapses)');
xlabel('f \times N^4');
ylabel('Convergence time');
leg = legend(labels, 'Location', 'northwest');
title(leg, 'Average connectivity');
hold off