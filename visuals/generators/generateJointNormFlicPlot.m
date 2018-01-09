% Loads data and creates a 2D scatter plot with logarithmic axes
clear variables
bkg = [ .9 .9 .9 ];

load results/cc_all.mat
X1 = (D_values.^regRes(2)) .* (N_values.^regRes(3));
Y1 = C_values;
Fval = D_values;

load results/fli01.mat
X2 = (D_values.^regRes(2)) .* (N_values.^regRes(3));
Y2 = C_values;
R = regress(Y2, X2);

Fval = [Fval; D_values];
X = [X1; X2];
Y = [Y1; Y2];

xMin = min(X);
xMax = max(X);
yMin = min(Y);
yMax = max(Y);
Xfit = logspace( log10(xMin), log10(xMax), 1000 );
Yfit1 = Xfit * regRes(1);
Yfit2 = Xfit * R;

xLab = ['f^{',num2str(regRes(2),3),'} \times N^{',num2str(regRes(3),3),'}'];
labs = { 'Tabula Rasa', 'Flickering', ...
    [ num2str( regRes(1),3),  '\times' , xLab ], ...
    [ num2str(R,3), ' \times ', xLab ] };

% Fig 1 corresponds to different colors densities
figure(10);
clf;
whitebg(bkg);
colormap parula;
scatter( X1, Y1, ...
    30, [.1 .4 .7], 'o', 'filled', 'MarkerFaceAlpha', 0.8 );
set( gca, 'xscale', 'log', 'yscale', 'log');
hold on;
scatter( X2, Y2, ...
    30, [.7 .4 .1], 'o', 'filled', 'MarkerFaceAlpha', 0.8 );
set( gca, 'xscale', 'log', 'yscale', 'log');
nFit = plot(Xfit,Yfit1,':k','LineWidth',2.5);
fFit = plot(Xfit,Yfit2,'-.k','LineWidth',2.5);
hold off;

legend( labs, 'location', 'northwest' );
axis( [xMin xMax yMin yMax] );
title('Convergence');
xlabel([ ...
    'd^{', num2str(regRes(2),3), '} \times ', ...
    'N^{', num2str(regRes(3),3), '}']);
ylabel('Convergence time');
set(gca,'FontSize',20)
