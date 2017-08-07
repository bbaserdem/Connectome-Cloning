% Use this function to generate two graphs, one showing a function with
% local minima that are not global, and one with equivelant minima
clear variables

% Number of minima to demonstrate
M = 5;
% Points per hump
T = 50;

% Function size
S = 3;

% Generate X values
Xval = 0:(pi/T):(M*2*pi);
Mind = 1 + T:(2*T):(2*M*T);

% Generate `true` function
fun = S * ( cos(Xval) + 1 ) / 2;
Tbas = zeros( size(Xval) );
TX = Xval( Mind );
TY = zeros( size(TX) );

Fbas = imgaussfilt( randn(size(Xval)), 30 );
FX = TX;
FY = fun(Mind) + Fbas(Mind);
[FYmin, FXmin] = min(FY);
FX(FXmin) = [];
FY(FXmin) = [];
FXmin = TX(FXmin);


figure(1);
subplot(1,2,1);
plot(Xval,fun+Tbas,Xval,Tbas,':','LineWidth',2);
title('Function with equivelant local minima')
hold on
scatter(TX,TY,200,...
    'MarkerFaceColor', [.3, .9, .3],...
    'MarkerEdgeColor', [.1, .5, .2],...
    'MarkerFaceAlpha', 0.5);
hold off
legend({'Function', 'Minima Surface', 'Global Minima'});
axis([Xval(1), Xval(end), -0.5, S+0.5]);

subplot(1,2,2)
plot(Xval,fun+Fbas,Xval,Fbas,':','LineWidth',2);
title('Function with non-equivelant local minima');
hold on
scatter(FX,FY,150,...
    'MarkerFaceColor', [.3, .4, .3],...
    'MarkerEdgeColor', [.1, .2, .2],...
    'MarkerFaceAlpha', 0.5);
scatter(FXmin,FYmin,200,...
    'MarkerFaceColor', [.3, .9, .3],...
    'MarkerEdgeColor', [.1, .5, .2],...
    'MarkerFaceAlpha', 0.5);
hold off
legend({'Function', 'Minima Surface', 'Local Minima','Global Minima'});
axis([Xval(1), Xval(end), -0.5, S+0.5]);
