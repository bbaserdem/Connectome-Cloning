function [ connectome ] = cc_genConnectome( N, d )
% CC_genConnectome: Generates a connectome.
%   N: Size of connectome
%   d: Density. Dictates algorithm to use.
%       - d>=N generates all-to-all network.
%       - If d >= 1, generates a scale-free 3rd power law scaling network
%           using the Barabasi-Albert algorithm. (d is rounded to integer)
%       - If 0<d<1 generates a random network (Erdos algorithm).
%       - d<=0 generates empty network.
%   connectome: The resulting connectome, given as a matrix of 1's and 0's.

if d >= N
    connectome = ones(N);
elseif d >= 1
    % Barabasi-Albert algorithm.
    % Initialize empty network with all-to-all core.
    connectome = zeros(N);
    d = round(d);
    connectome( 1:d , 1:d ) = 1;
    new= zeros(1, d);
    % Add nodes with preferential attachment
    for a = d+1:N
        new(:) = 0;
        pro = [ cumsum( sum( ...
            connectome( (1:(a-1)), (1:(a-1)) ) + (...
            connectome( (1:(a-1)), (1:(a-1)) )' ) ...
            ) ) / (2*d*(a-1) + 1), 1];
        % find first added connection
        new(1) = find( pro > rand, 1);
        % add the other connections
        for b = 2:d
            while   (  new(b)==0) || ...
                    ( (new(b)~=a) && ( sum(new == new(b)) > 2 ) ) || ...
                    ( (new(b)==a) && ( sum(new == new(b)) > 1 ) )
                % line1; pick a node if it will be the first choice
                % line2; pick a node if there are already 2 of the same con
                % line3; pick a node if self-con already selected before
                new(b) = find( pro > rand, 1);
            end
        end
        % add connections to connectome(random direction)
        for b = 1:d
            % If a twice repeated index, add the second direction
            if     connectome( new(b),      a) == 1
                connectome(         a, new(b)) = 1;
            elseif connectome(      a, new(b)) == 1
                connectome(    new(b),      a) = 1;
            % If empty, randomly assign direction
            elseif rand < 0.5
                connectome(    new(b),      a) = 1;
            else
                connectome(         a, new(b)) = 1;
            end
        end
    end
elseif d > 0
    % Erdos random initialization
    connectome = rand(N) < d;
else
    connectome = zeros(N);
end

end