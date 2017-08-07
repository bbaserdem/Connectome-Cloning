function [ C, f ] = init_conn( N, method, p )
%INIT_CONN Function to generate networks of different type
%   Network types programmed;
%   1 > 'Erdos':    (Default) Network by random connectivity. 
%   2 > 'Barabasi': Barabasi-Albert mdoel for scale free
%   C:      Returns sparse connectivity graph
%   e2e:    Returns index list of connectivity
%   fn4:    Returns usededges*possibleedges
% Defaults to density 0.2 for erdos, and m0=10, m=1 for Barabasi

% Default to erdos
if ~exist('method', 'var')
    method = 'Erdos';
end

if strcmp(method,    'Erdos')
    if ~exist('p', 'var')
        p = 0.2;
    elseif p > 1
        error('Invalid density for Erdos graph');
    end
    C = rand(N) < p;
    f = sum( C(:) );
elseif strcmp(method, 'Barabasi')
    if ~exist('p', 'var')
        p = 2;
    end
    if p > N
        error('Invalid seed size, %d(N) > %d(m0)\n', N, p(1));
    end
    % Seed network
    C = ones(p);
    C(N,N) = 0;
    
    % Grow nodes with preferantial attachment
    for n = (p+1):N
        % generate cumprob
        adds = zeros(1,p);
        prob = sum( C(1:(n-1),1:(n-1)) + C(1:(n-1),1:(n-1))' );
        prob = [cumsum( prob / ( 2*p*(n-1) + 1 ) ), 1];
        % keep a list of connection to add to that node
        adds(1) = find( prob > rand, 1);    % First step
        for m = 2:p
            % Add each node randomly (can add twice for a multidir node)
            while   (  adds(m)==0) || ...
                    ( (adds(m)~=n) && ( sum(adds == adds(m)) > 2 ) ) || ...
                    ( (adds(m)==n) && ( sum(adds == adds(m)) > 1 ) )
                % do uniquely
                adds(m) = find( prob > rand, 1);
            end
        end
        % add nodes with random direction
        for m = 1:p
            % If a twice repeated index, add the second direction
            if     C(adds(m),n) == 1
                C(n,adds(m)) = 1;
            elseif C(n,adds(m)) == 1
                C(adds(m),n) = 1;
            % If empty, randomly assign direction
            elseif rand < 0.5
                C(adds(m),n) = 1;
            else
                C(n,adds(m)) = 1;
            end
        end
    end
    f = p * N;
else
    error('Invalid generation option');
end

end

