% Loads data and creates a 2D scatter plot with logarithmic axes
clear variables
load results/all.mat

% Turn on to generate also a figure with exponential fit
exponentialFit = true;
bkg = [ .9 .9 .9 ];

X = (D_values.^regRes(2)) .* (N_values.^regRes(3));
Y = C_values;
xMin = min(X);
xMax = max(X);
yMin = min(Y);
yMax = max(Y);
Xfit = logspace( log10(xMin), log10(xMax), 1000 );
Yfit = Xfit * regRes(1);

% Fig 1 corresponds to different colors densities
figure(1);
clf;
whitebg(bkg);
colormap parula;
scatter( X, Y, ...
    40, log10( D_values ), 'o', 'filled', 'MarkerFaceAlpha', 0.8 );
set( gca, 'xscale', 'log', 'yscale', 'log');
hold on;
plot(Xfit,Yfit,':k','LineWidth',1.5);
hold off;
axis( [xMin xMax yMin yMax] );
title('Convergence');
xlabel([ ...
    'd^{', num2str(regRes(2),3), '} \times ', ...
    'N^{', num2str(regRes(3),3), '}']);
ylabel('Convergence time');
barTemp = colorbar;
title(barTemp,'(log_{10}) Network Density');
set(gca,'FontSize',30)

% Fig 2 corresponds to different network sizes
figure(2);
clf;
whitebg(bkg);
colormap jet;
scatter( X, Y, ...
    40, log10( N_values ), 'o', 'filled', 'MarkerFaceAlpha', 0.8 );
set( gca, 'xscale', 'log', 'yscale', 'log');
hold on;
plot(Xfit,Yfit,':k','LineWidth',1.5);
hold off;
axis( [xMin xMax yMin yMax] );
title('Convergence');
xlabel([ ...
    'd^{', num2str(regRes(2),3), '} \times ', ...
    'N^{', num2str(regRes(3),3), '}']);
ylabel('Convergence time');
barTemp = colorbar;
title(barTemp,'(log_{10}) Network Size');
set(gca,'FontSize',30)

% Fig 3 & 4 will correspond to exponential fit (will do later)

% run exponential regression ( C = X1 * exp(D*X2) * N^X3
regExp = regress( ...
    log(C_values), ...
    [ ones(size(C_values)), D_values, log(N_values) ] );
regExp(1) = exp(regExp(1));
E = exp(D_values*regExp(2)) .* (N_values.^regRes(3));
eMin = min(E);
eMax = max(E);
Efit = logspace( log10(eMin), log10(eMax), 1000 );
YEfit = Efit * regExp(1);

% Fig 3 corresponds to different colors densities
figure(3);
clf;
whitebg(bkg);
colormap parula;
scatter( E, Y, ...
    40, log10( D_values ), 'o', 'filled', 'MarkerFaceAlpha', 0.8 );
set( gca, 'xscale', 'log', 'yscale', 'log');
hold on;
plot(Efit,YEfit,':k','LineWidth',1.5);
hold off;
axis( [eMin eMax yMin yMax] );
title('Convergence');
xlabel([ ...
    'e^{', num2str(regExp(2),3), '\times d}\times', ...
    'N^{', num2str(regExp(3),3), '}']);
ylabel('Convergence time');
barTemp = colorbar;
title(barTemp,'(log_{10}) Network Density');
set(gca,'FontSize',30)

% Fig 4 corresponds to different network sizes
figure(4);
clf;
whitebg(bkg);
colormap jet;
scatter( E, Y, ...
    40, log10( N_values ), 'o', 'filled', 'MarkerFaceAlpha', 0.8 );
set( gca, 'xscale', 'log', 'yscale', 'log');
hold on;
plot(Efit,YEfit,':k','LineWidth',1.5);
hold off;
axis( [eMin eMax yMin yMax] );
title('Convergence');
xlabel([ ...
    'e^{', num2str(regExp(2),3), '\times d}\times', ...
    'N^{', num2str(regExp(3),3), '}']);
ylabel('Convergence time');
barTemp = colorbar;
title(barTemp,'(log_{10}) Network Size');
set(gca,'FontSize',30)




























