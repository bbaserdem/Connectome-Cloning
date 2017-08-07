N = ceil( logspace(1,1.5,5) );
D = 0.1:0.1:0.9;
repeat = 10;

density_con       = cell( length(N), length(D) );
convergence       = cell( length(N), length(D) );
[density_con{:}]  = deal( zeros(1,     repeat) );
[convergence{:}]  = deal( zeros(1,     repeat) );

fprintf('Doing jittering network\n');
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
            [convergence{n,d}(r), density_con{n,d}(r)] = ...
                conclone_wrapper( nt, dt );
        end
        toc;
    end
end

% Determine coefficients
X = [];
Y = [];
for d = 1:length(D)
    for n = 1:length(N)
        index = convergence{n,d}~=0;
%         X = [ X; ...
%             ones(sum(index),1), ...
%             density_con{n,d}(index)', ...
%             log( N(n) ) * ones( sum(index),1 )   ];
%         Y = [ Y; log( convergence{n,d}(index) )' ];
        X = [ X; ...
            ones(sum(index),1), density_con{n,d}(index)' ];
        Y = [ Y; log( convergence{n,d}(index) )' - 4*log(N(n)) ];
    end
end
B = regress(Y,X);

% Start plotting
colors = hsv2rgb( [...
    (1:length(D))'/length(D), 0.9*ones(length(D),2) ] );
labels = cell(1, length(D) );

figure(1);
hold on;
for d = 1:length(D)
    labels{d} = num2str(D(d));
    x = [];
    y = [];
    for n = 1:length( N )
        index = convergence{n,d}~=0;
%         x = [ x, density_con{n,d}(index).^( B(2) ) * (N(n)^B(3)) ];
%         x = [ x, exp( B(2) * density_con{n,d}(index) ) * (N(n)^B(3))];
        x = [ x, exp( B(2) * density_con{n,d}(index) ) * (N(n)^4)];
        y = [ y, convergence{n,d}(index) ];
    end
    scatter(...
        x,...
        y,...
        20,...
        colors(d,:), 'filled' ...
        );
end
set( gca, 'xscale', 'log', 'yscale', 'log');
x = xlim;
x = logspace( log10(x(1)), log10(x(2)), 1000 );
y = x * exp(B(1));
plot(x,y,':k','LineWidth',1.5);
title('Convergence (Flickering synapses)');
% xlabel( ['e^{', num2str(B(2),3), '\timesf} \times N^{',...
%     num2str(B(3),3), '}'] );
xlabel( ['e^{', num2str(B(2),3), '\timesf} \times N^4'] );
ylabel('Convergence time');
leg = legend(labels, 'Location', 'northwest');
title(leg, 'Average connectivity');
hold off

save jitter_results