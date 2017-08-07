% Generate a gif of a sample run
clear variables
N = 1000;
D = 0.05;
fname = ['visuals/video/run', num2str(N), '.avi'];

% Want a better background;         [ 0.1, 0.1, 0.1 ]   gray
% Original connectome will be       [ 0.0, 0.8, 0.1 ]   green
% Reconstruction connectome will be [ 0.8, 0.0, 0.1 ]   magenta
% Correct will be then be           [ 0.9, 0.9, 0.3 ]   yellow
col = [
    .1 .1 .1
    .1 .4 .7
    .7 .4 .1
    .6 1. .6
    ];

% Run simulation
tic;
fprintf('Running simulation\n');
[ ~, oriConn, map, proc, stepSize] = cc_video( N, D );
F = size(proc,3);
toc;
fprintf('Simulation done!\n');

% Check if the thing is correct
if graphisomorphism( sparse(proc(:,:,end)), sparse(oriConn) )
    fprintf('The graphs are isomorphic\n');
end

% Open and setup videowriter
videoFile = VideoWriter( fname, 'Uncompressed AVI' );
% Shoot for 10 seconds
videoFile.FrameRate = max(1,floor(F/10));
open( videoFile );

% Determine if image needs resizing if too small (to >500x500)
% Assuming presentation will be 720p
scale = 1;
% if N <= ( 500 / 2 )
%     scale = floor( 500 / N );
%     fprintf('Will scale image by factor %d to suit 720p.\n', scale);
%     fprintf('Resulting video size will be %dx%d.\n', N*scale, N*scale);
% elseif N > 500
%     fprintf('Resulting video is too large for 720p presentation;\n');
%     fprintf('Will downscale to 500x500. Quality will decrease.');
%     scale = 500 / N;
% else
%     fprintf('Resulting video size will be %dx%d.\n', N, N);
% end

% write each frame
frames = zeros(N,N,3,F);
for f = 1:F
    for rgb = 1:3
        frames(:,:,rgb,f) = ...
            col(1,rgb) * ~( proc(:,:,f) | oriConn ) + ...   % Background
            col(2,rgb) * ( oriConn & ~proc(:,:,f) ) + ...   % Missing
            col(3,rgb) * ( ~oriConn & proc(:,:,f) ) + ...   % False
            col(4,rgb) * ( oriConn & proc(:,:,f)  );        % True
    end
end
if scale > 1
    frames = imresize( frames, scale, 'nearest' );
elseif scale < 1
    frames = imresize( frames, scale );
end
writeVideo( videoFile, frames );

% Close video file
close( videoFile );
fprintf('Done with video!\n');