% Loads data and creates a 2D scatter plot with logarithmic axes
clear variables
load results/fli01.mat

% Turn on to generate also a figure with exponential fit
% exponentialFit = true;
bkg = [ .9 .9 .9 ];

% Generate regression
[R, Rconf, ~, ~, Rstats] = regress( ...
    log( C_values ), ...
    [ ones(size(C_values)), log(D_values), log(N_values) ] );
R(1) = exp( R(1) );
Rconf(1,:) = exp( Rconf(1,:) );

X = (D_values .^ R(2) ) .* ( N_values .^ R(3) );
Y = C_values;
xMin = min(X);
xMax = max(X);
yMin = min(Y);
yMax = max(Y);
Xfit = logspace( log10(xMin), log10(xMax), 1000 );
Yfit = Xfit * R(1);

figure(5);
clf;
whitebg(bkg);
colormap parula;

% Plot densities.
subplot(1,2,1);
scatter( X, Y, ...
    40, log10( D_values ), 'o', 'filled', 'MarkerFaceAlpha', 0.8 );
set( gca, 'xscale', 'log', 'yscale', 'log');
hold on;
store = plot(Xfit,Yfit,':k','LineWidth',1.5);
legend(store,[...
    'C = ', num2str(R(1),3), ...
    ' d^{', num2str(R(2),3), ...
    ' }N^{', num2str(R(3),3), '}' ], 'location', 'southeast');
hold off;
axis( [xMin xMax yMin yMax] );
title('Convergence (Flickering Synapses)');
xlabel([ ...
    'd^{', num2str(R(2),3), '} \times ', ...
    'N^{', num2str(R(3),3), '}']);
ylabel('Convergence time');
barTemp = colorbar;
title(barTemp,'Network Density_{log_{10}}');
set(gca,'FontSize',15);

% Show failure amount
subplot(1,2,2);
imagesc(failure);
barTemp = colorbar;
title(barTemp,'Failure rate');
title( 'Failing simulations' );
xlabel( 'Average network connection density' );
ylabel( 'Number of neurons in network' );
set(gca, 'YTickLabel', N, 'XTickLabel', round(D,2), 'FontSize', 15 );
