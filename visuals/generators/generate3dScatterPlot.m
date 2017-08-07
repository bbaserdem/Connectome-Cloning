% Loads data and creates a 3D scatter plot with logarithmic axes
clear variables
load results/all.mat
% Contains N_values, C_values, D_values, regressionResults

% Grab some info
n_min = min( N_values );
n_max = max( N_values );
d_min = min( D_values );
d_max = max( D_values );
c_min = min( C_values );
c_max = max( C_values );

% Colors for differentiating values are in an hsv range.
bkgr = [.9, .9, .9];
color_style = 'dens';
hueSat = [1, 0.5];

% Color matrix for scatterplot
if      strcmp( color_style, 'dens' )
    cols = ( log(D_values) - log(d_min) ) / (log(d_max)-log(d_min));
elseif  strcmp( color_style, 'size' )
    cols = ( log(N_values) - log(n_min) ) / (log(n_max)-log(n_min));
elseif  strcmp( color_style, 'conv' )
    cols = ( log(C_values) - log(c_min) ) / (log(c_max)-log(c_min));
else
    error('Invalid color method');
end
cols = hsv2rgb( [cols, repmat(hueSat, size(cols,1), 1) ] );

% Generate the scatter plot
whitebg(bkgr);
figure(1);
clf;
scatter3( D_values, N_values, C_values, ...
    50, cols, 'filled', 'MarkerFaceAlpha', 0.5 );
title('Convergence time vs. Network Parameters');
xlabel('Network density');
ylabel('Network size');
zlabel('Convergence time');
set( gca, ...
    'xscale', 'log', 'yscale', 'log', 'zscale', 'log', ... 
    'XMinorTick', 'off', 'YMinorTick', 'off', 'ZMinorTick', 'off' );

% Draw the surface on the plot area
[X,Y] = meshgrid(...
    logspace( log10(d_min), log10(d_max), 20), ...
    logspace( log10(n_min), log10(n_max), 20) );
Z = regRes(1) .* ...
    ( X .^ regRes(2) ) .* ...
    ( Y .^ regRes(3) );
colormap summer;
hold on;
surf( X, Y, Z, log10(Z), ...
    'FaceAlpha', 0.4, ...
    'EdgeColor', 'none');
hold off;
set(gca,'FontSize',24)