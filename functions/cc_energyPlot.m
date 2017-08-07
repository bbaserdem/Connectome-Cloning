function [ Xval, Yval ] = cc_energyPlot( N, D, recordAmount )
%CONNCLONE Simulation of a connectome cloning run (added energy recording)
%   Build on top of the function as of 15-07-2017
%   Uses results from the previous run to estimat runtime

% Inputs
%   N: ( 100 ) input connectome size
%   D: ( 0.1 )input connectome density. >=1 switches to scale free k=3
%       network, which then becomes half the average degree per node
%   E: (  10 ) Cost parameter in hamiltonian
%   H: ([0,0]) Temperature, scaling from H(1) to H(2)
%
% Outputs are video specific
%   D: density of connectivity
%   C: Original connectome
%   conid: Mapping between original and reconstructed
%   VID: 3D logical matrix, with the connectome of reconstruction
%   stepSize: Number of simulation steps between each step
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
BAT = 500000;
E = 10;
H = [0, 0];

%-----PARAMETERS-----% END






%% %-----INITIALIZATION-----%

% For video 
cEstimate = round( exp(2.5199) * ( D ^ 1.6309 ) * ( N ^ 3.3652 ) );
% Aim to record ~10000 recordings
stepSize = ceil( cEstimate / recordAmount );
% Account for overhead
Xval = (1:stepSize:2*cEstimate)';
Yval = zeros(size(Xval));
lEstimate = length(Xval);

% Connectome
C = cc_genConnectome(N, D);      % Generate connectome
B = sum( C(:) );
D = B / (N^2);
% Time steps
T = round( 5 * D * (N^4) );
% Move probabilities ( jump (1-D), swap(D-1/N), flip(1/N) ) (cumsum)
R = [ 1-D, max(1-D,1-1/N) ];

% Function to evaluate transition probability
if all( H == 0 )
    tran = @(t,de) 1*(de<0);
else
    tran = @(t,de) exp( -1 * de / (H(1) + ((t-2)/(T-2))*(H(2)-H(1)) ) );
end

% Parent connectome barcode-list generation
[K(:,1), K(:,2), ~] = find( C );
% (Barcode-pair,synapse) list, initialized randomly
S = ceil( N * rand(B, 2) );

% Barcode matrix, (barcode,cell) and populate it from S
M = zeros(N);
J = zeros(N,1);
A = zeros(N,B);         % BP list in each cell
Q = zeros(N,B);         % Positions in A
for b = 1:B
    M( K(b,1), S(b,1)) = M( K(b,1), S(b,1)) + 1;
    M( K(b,2), S(b,2)) = M( K(b,2), S(b,2)) + 1;
    J( S(b,1) ) = J(S(b,1)) + 1;
    A( S(b,1), J(S(b,1)) ) = b;
    Q( S(b,1),         b ) = J(S(b,1));
    % Just once if on a self synapse
    if S(b,1) ~= S(b,2)
        J( S(b,2) ) = J( S(b,2) ) + 1;
        A( S(b,2), J(S(b,2)) ) = b;
        Q( S(b,2),         b ) = J( S(b,2) );
    end
end

% Number of barcodes per cell
F = sum(M);

% Store cells with more than one barcode to speed simulation
Y = zeros(1,N);
X = zeros(1,N);
Z = 0;
for n = 1:N
    if J(n) > 1
        Z = Z + 1;
        Y(Z) = n;
        X(n) = Z;
    end
end

% Energy function
L = ...
    sum( -(1+E)*sum(M.^2) + E*(sum(M).^2) ) + ...   % Current E contrib.
    sum( sum(C+C').^2 );                % Normalize GS

%-----INITIALIZATION-----% END






%% %-----SIMULATION-----%
step = 0;
for t = 1:T
    
    %-----Energy stuff!-----%
    if mod(t-1, stepSize) == 0
        step = step + 1;
        Yval(step) = L;
    end
    
    k = mod( t-1, BAT ) + 1;
    if k == 1
        r = rand(5,BAT);
    end
    
    %-----ACTION-----% BEGIN
    if r(1,k) < R(1)
        %-----JUMP-----% BEGIN
        b = ceil( B * r(2,k) );   % Select a barcode
        p = ceil( 2 * r(3,k) );   % Select an end
        c  = S(b,3-p);          % Base cell
        c1 = S(b,p);            % Changing cell is S(b,p)
        c2 = ceil( N * r(4,k) );	% Changed cell is a random cell not c1
        if c1 == c2
            continue
        end
        m = K(b,p);             % The barcode that moves is K(b,p)
        
        % Change in energy
        dE = 2 * ( E*(1+F(c2)-F(c1)) - (1+E)*(1+M(m,c2)-M(m,c1)) );
        
        % Accept clause
        if r(5,k) < tran(t,dE)
            % If jumping from a self synapse, dont do anything
            if c ~= c1
                % Remove this barcode-pair from c1
                if J(c1) ~= 1
                    A(c1, Q(c1, b) ) = A(c1, J(c1) );
                    Q(c1, A(c1, J(c1)) ) = Q(c1, b);
                    Q(c1, b) = 0;
                    A(c1, J(c1)) = 0;
                else
                    A(c1, 1) = 0;
                    Q(c1, b) = 0;
                end 
                J(c1) = J(c1) - 1;
                % If cell is left with only 1 barcode-pair;
                if J(c1) == 1 % remove cell from Y
                    Y( X(c1) ) = Y( Z );
                    X( Y(Z) ) = X(c1);
                    Y( Z ) = 0;
                    X( c1 ) = 0;
                    Z = Z - 1;
                end
            end
            % If jumping to a self synapse, dont do anything
            if c ~= c2
                % Add this barcode-pair to c2
                J(c2) = J(c2) + 1;
                A(c2, J(c2) ) = b;
                Q(c2,     b ) = J(c2);
                % If a cell now has 2 barcode pairs, add the cell to Y
                if J(c2) == 2
                    Z = Z + 1;
                    Y( Z ) = c2;
                    X( c2 ) = Z;
                end
            end
            F(c1) = F(c1) - 1;
            F(c2) = F(c2) + 1;
            M(m,c1) = M(m,c1) - 1;
            M(m,c2) = M(m,c2) + 1;
            S(b,p) = c2;
            L = L + dE;
        end
        %-----JUMP-----% END
    elseif r(1,k) < R(2)
        %-----SWAP-----% BEGIN
        % Pick cell with 2 barcode-pairs
        c = Y( ceil( Z * r(2,k) ) );   
        % Pick b1 & b2 randomly among barcodes inside the cell
        b1 = A(c, ceil( J(c) * r(3,k) ) );
        b2 = A(c, ceil( J(c) * r(4,k) ) );
        if b1 == b2
            continue
        end
        
        % Get p1 and p2 ends
        if     S(b1, 1) ~= c
            p1 = 1;
        elseif S(b1, 2) ~= c
            p1 = 2;
        else % assign randomly if self synapse
            p1 = ceil( 2 * rand );
        end
        if     S(b2, 1) ~= c
            p2 = 1;
        elseif S(b2, 2) ~= c
            p2 = 2;
        else % assign randomly if self synapse
            p2 = ceil( 2 * rand );
        end
        % Get the cells and barcodes
        c1 = S(b1,p1);
        c2 = S(b2,p2);
        m1 = K(b1,p1);
        m2 = K(b2,p2);
        % Energy difference
        if (c1==c2)
            continue
        elseif (m1==m2)
            dE = 0;
        else
            dE = -2 * (1+E) * ( 2 - (...]
                M(m1,c1) + M(m2,c2) ) + ...
                M(m1,c2) + M(m2,c1) );
        end
        % Accept clause
        if r(5,k) < tran(t,dE)
            % If c-c1 a self synapse, dont do anything
            if c ~= c1
                % Remove barcode b1 from, add b2 to c1
                A(c1, Q(c1, b1) ) = b2;
                Q(c1, b2) = Q(c1, b1);
                Q(c1, b1) = 0;
            end
            % If c-c2 a self synapse, dont do anything
            if c ~= c2
                % Remove barcode b2 from, add b1 to c1
                A(c2, Q(c2, b2) ) = b1;
                Q(c2, b1) = Q(c2, b2);
                Q(c2, b2) = 0;
            end
            M(m1,c1) = M(m1,c1) - 1;
            M(m2,c1) = M(m2,c1) + 1;
            M(m2,c2) = M(m2,c2) - 1;
            M(m1,c2) = M(m1,c2) + 1;
            S(b1,p1) = c2;
            S(b2,p2) = c1;
            L = L + dE;
        end
        %-----SWAP-----% END
    else
        %-----FLIP-----% BEGIN
        b = ceil( B * r(2,k) );           % Randomly choose a barcode-pair
        c1 = S(b,1);                    % c1 is the pre-synaptic end
        c2 = S(b,2);                    % c2 is the postsynaptic end
        m1 = K(b,1);                    % m1 is the pre-synaptic barcode
        m2 = K(b,2);                    % m2 is the postsynaptic barcode
        % Energy difference
        if c1 == c2 
            continue
        elseif m1 == m2
            dE = 0;
        else
            dE = -2 * (1+E) * ( 2 - (...]
                M(m1,c1) + M(m2,c2) ) + ...
                M(m1,c2) + M(m2,c1) );
        end
        % Accept clause
        if r(5,k) < tran(t,dE)
            M(m1,c1) = M(m1,c1) - 1;
            M(m2,c1) = M(m2,c1) + 1;
            M(m2,c2) = M(m2,c2) - 1;
            M(m1,c2) = M(m1,c2) + 1;
            S(b,1) = c2;
            S(b,2) = c1;
            L = L + dE;
        end
        %-----FLIP-----% END
    end
    %-----ACTION-----% END
    
    
    %-----TRUNCATE-----% BEGIN
    if L == 0
        %-----Energy stuff!
        step = step+1;
        Xval(step) = t;
        Xval(step+1:end) = 0;
        if step > lEstimate
            Xval = Xval(1:lEstimate);
            Yval = Yval(1:lEstimate);
        end
        
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
    Xval = []; Yval = [];    % Fail if not OBOC
else
    V = sparse( conid(S(:,1)), conid(S(:,2)), 1, N, N);
    if ~all( C(:)==V(:) )
        Xval = []; Yval = [];
    end
end
%-----RECONSTRUCTION-----% END

end

