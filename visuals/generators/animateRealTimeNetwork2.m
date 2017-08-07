% This script will animate a real time network.
% Uses function cc_realTimeConnectome
% This network is drawn a differently than the last one

clear variables

% Initial connectome
N = 15;
C = zeros(N);
S = [       % First module (pentagram)
    1   4
    2   5
    3   1
    3   5
    4   2 ];
S = [ S     % Second module (pentagon)
    6   6
    6   7
    6   10
    7   8
    9   8
    10  9 ];
S = [ S     % Third module (cross stitch)
    11  13
    12  11
    13  11
    13  15
    14  14
    15  13
    15  14 ];
S = [ S     % Intermodule connections
    2   14
    7   12
    8   11
    8   1
    10  5
    11  8
    12  7
    14  3
    14  2 ];
    
for k = 1:size(S,1)
    C(S(k,1),S(k,2)) = 1;
end
graphONN = digraph(C);
names = cell(N,1);
for k = 1:N
    names{k} = ['Cell ', char(64+k)];
end
graphONN.Nodes.Name = names;
R = 5;
n = 5; % submodule size
c = 3; % submodule number

X = repmat( R * sin( (1:n) * 2 * pi / n ), 1, c) + ...
    kron( -2*R * sin( (1:c) * 2 * pi / c ), ones(1,n) );
Y = repmat( R * cos( (1:n) * 2 * pi / n ), 1, c) + ...
    kron( -2*R * cos( (1:c) * 2 * pi / c ), ones(1,n) );

% Draw connectome
f = figure(1);
set(gcf,'position', [100 100 1400 700 ] );
subplot(1,2,1);
whitebg( [.8 .8 .8] );
plot(graphONN, 'XData', X, 'YData', Y, ...
    'LineWidth', 2, 'ArrowSize', 15, 'EdgeColor', [.4 .6 .7] );
title('Original connectome')
axis( 3.5 * [-R R -R R] );
axis off;

% Run reconstruction
[aniCon, map] = cc_realTimeConnectome( C );
T = size(aniCon, 3);
dig = floor( log10(T) ) + 1;
% Rearrange positions
X = X(map);
Y = Y(map);
for k = 1:N
    names{k} = ['Cell ', num2str(k) ];
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
    graphConn.Nodes.Name = names;
    plot(graphConn, 'XData', X, 'YData', Y, ...
        'LineWidth', 3, 'ArrowSize', 20, 'EdgeColor', [.7 .4 .1] );
    title( sprintf( [ 'Step: %', num2str(dig), 'd(/%d)'], t, T ) );
    axis( 3.5 * [-R R -R R] );
    axis off;
    movie(p) = getframe(gcf);
end

% Save video
vid = VideoWriter( 'visuals/video/realtime15', 'Uncompressed AVI' );
% Compress to ~ 10 seconds
vid.FrameRate = floor(p/10);
open(vid);
writeVideo(vid, movie);
close(vid);








































