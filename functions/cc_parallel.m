function [ O, D ] = cc_parallel( N, D, E )
%CC_PAR Parallel version of the cloning script
%   This version tries to jump every barcode pair simultaneously on one it.

% Inputs
%   N: ( 100 ) input connectome size
%   D: ( 0.1 )input connectome density. >=1 switches to scale free k=3
%       network, which then becomes half the average degree per node
%   E: (  10 ) Cost parameter in hamiltonian
%   H: ([0,0]) Temperature, scaling from H(1) to H(2)
%
% Outputs
%   O: 0 if reconstruction was unsuccessful, convergence time if
%       successful.
%   B: number of barcode pairs in the connectome
%
% Intermediate varaibles
%   BAT: Used to batch generate random numbers
%   C: The original connectome
%   V: Reconstructed connectome
%   T: Max simulation runtime
%   R: Move probabilities
%   K: List of barcodes on barcode-pair(bp) ends (pre,pos) 
%   S: List that contains bp locations. (pre,pos)
%   M: Table keeping track of barcodes per cell (barc,cell)
%   F: Total barcodes per cell (cell)
%   A: List of bp ids in each cell (cell,...)
%   J: Total # of barcode pairs per cell (cell)
%   Q: Location of bp in A (cell,bp-id)
%   Y: List of cells with >1 bp in them (...)
%   Z: Total # of cells with >1 bp in them
%   X: Location of cells in A (cell)
%   c, c1, c2: used to mark selected cells
%   b, b1, b2: used to mark selected barcode-pair ids
%   p, p1, P2: used to mark pre-pos ends of selected bps
%   m, m1, m2: used to mark barcodes





%% %-----PARAMETERS-----%
BAT = 1000;
if ~exist('N', 'var')
    N = 100;
end
if ~exist('D', 'var')
    D = 0.1;
end
if ~exist('E', 'var')
    E = 10;
end

%-----PARAMETERS-----% END






%% %-----INITIALIZATION-----%

% Connectome
C = cc_genConnectome(N, D);      % Generate connectome
B = sum( C(:) );
D = B / (N^2);
% Time steps
T = 10*cc_averageRuntime(N,D);

% Parent connectome barcode-list generation
[K(:,1), K(:,2), ~] = find( C );
% (Barcode-pair,synapse) list, initialized randomly
S = ceil( N * rand(B, 2) );

% Barcode matrix, (barcode,cell) and populate it from S
M = zeros(N); 
for b = 1:B
    M( K(b,1), S(b,1)) = M( K(b,1), S(b,1)) + 1;
    M( K(b,2), S(b,2)) = M( K(b,2), S(b,2)) + 1;
end
% Number of barcodes per cell
F = sum(M);

% Ground state
GS = -1 * sum( sum(C+C').^2 );

m = zeros(1,B);
c = zeros(1,B);
d = zeros(1,B);
V = 0;
%-----INITIALIZATION-----% END






%% %-----SIMULATION-----%

for t = 1:T
    k = mod( t-1, BAT ) + 1;
    if k == 1
        % Preselect barcode ends and stuff
        R = rand(BAT,2*B);
    end
    
    %-----JUMPS-----% BEGIN
    
    % new cell candidates
    n = ceil( N * R(k,1:B) );
    % new positions
    p = ceil( 2 * R(k,((1:B)+B)) );
    % Calculate dE
    for b = 1:B
        % Cell to move is n(b)
        % Synaptic end to move is p(b)
        % Get the moved barcode
        m(b) = K(b, p(b) );
        % Get the origin cell
        c(b) = S(b, p(b) );
        % Calculate e difference
        if c(b) == n(b)
            d(b) = 0;
        else
            d(b) = 2 * ( E * ( 1 + F(n(b)) - F(c(b)) ) ...
                - (1+E) * ( 1 + M(m(b),n(b)) - M(m(b),c(b)) ) );
        end
    end
    % Accept moves if wanted
    for b = 1:B
        if d(b) <= 0
            S(b, p(b)) = n(b);
            M( m(b), c(b) ) = M( m(b), c(b) ) - 1;
            F( c(b) ) = F( c(b) ) - 1;
            M( m(b), n(b) ) = M( m(b), n(b) ) + 1;
            F( n(b) ) = F( n(b) ) + 1;
        end
    end
    %-----JUMPS-----% END
    
    % Energy
    L = sum( -(1+E)*sum(M.^2) + E*(sum(M).^2) ) -GS;
    
    %-----TRUNCATE-----% BEGIN
    if L == 0
        T = t;
        break
    end
    %-----TRUNCATE-----% END
end

%-----SIMULATION-----% END



%% %-----RECONSTRUCTION-----%

% Assign each cell a number that indicates their barcode ID
% Identify cells by most prominant barcode
conid = cc_isOBOC( M );
% check for repeating barcode-id's
if conid == false
    O = 0;    % Fail if not OBOC
else
    V = sparse( conid(S(:,1)), conid(S(:,2)), 1, N, N);
    if all( C(:)==V(:) )
        O = T;
    else
        O = 0;
    end
end
%-----RECONSTRUCTION-----% END

end

