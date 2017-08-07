% This script will animate a real time network.
% Uses function cc_realTimeConnectome

clear variables

% Initial connectome
N = 10;
C1 = zeros(N);
S1 = [
    1 1
    1, 6
    2, 10
    3, 4
    4, 8
    5, 2
    6, 1
    7, 5
    7, 7
    8, 9
    9, 3
    9, 8
    10, 7 ];
for k = 1:size(S1,1)
    C1(S1(k,1),S1(k,2)) = 1;
end
graphOri1 = digraph(C1);
names1 = cell(N,1);
for k = 1:N
    names1{k} = ['Cell ', char(64+k)];
end
graphOri1.Nodes.Name = names1;
% Generate circular coordinates for this graph, like a clock
R1 = 5;
X1 = (R1-1) * sin( (1:N) * 2 * pi / N );
Y1 = (R1-1) * cos( (1:N) * 2 * pi / N );
% Generate labels

% Draw connectome
f = figure(1);
set(gcf,'position', [100 100 1400 700 ] );
subplot(1,2,1);
whitebg( [.8 .8 .8] );
plot(graphOri1, 'XData', X1, 'YData', Y1, ...
    'LineWidth', 3, 'ArrowSize', 20, 'EdgeColor', [.4 .6 .7] );
title('Original connectome')
axis( [-R1 R1 -R1 R1] );
axis off;

% Run reconstruction
[aniCon, map] = cc_realTimeConnectome( C1 );
T = size(aniCon, 3);
dig = floor( log10(T) ) + 1;
% Rearrange positions
X1 = X1(map);
Y1 = Y1(map);
for k = 1:N
    names1{k} = ['Cell ', num2str(k) ];
end

% Generate movie
subplot(1,2,2);
p = 1;
for t = 1:T
    if t == 1
        old = aniCon(:,:,t);
    else
        if all( all( old == aniCon(:,:,t) ) )
            continue    % Go to next stage if no change occured
        else
            old = aniCon(:,:,t);
            p = p+1;
        end
    end
    graphConn = digraph( aniCon(:,:,t) );
    graphConn.Nodes.Name = names1;
    plot(graphConn, 'XData', X1, 'YData', Y1, ...
        'LineWidth', 3, 'ArrowSize', 20, 'EdgeColor', [.7 .4 .1] );
    title( sprintf( [ 'Step: %', num2str(dig), 'd(/%d)'], t, T ) );
    axis( [-R1 R1 -R1 R1] );
    axis off;
    movie(p) = getframe(gcf);
end

% Save video
vid = VideoWriter( 'visuals/video/realtime10', 'Uncompressed AVI' );
% Compress to ~ 10 seconds
vid.FrameRate = floor(p/10);
open(vid);
writeVideo(vid, movie);
close(vid);










































