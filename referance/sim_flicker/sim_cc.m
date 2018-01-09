function [ success, H ] = sim_cc( C,params )
%SIM_JO_G2 Simulation with quadratic hamiltonian (G=2)
% Inputs
%   C: input connectome
% Outputs
%   success: 0 if reconstruction was unsuccessful, convergence time if
%       successful.
%   conn: cell array, with conn{1} original and conn{2} reconstructed
%       connectome.
%   hamil: values of hamiltonian across time steps. Defaults to 10000
%       across projected runtime if larger.





%% %-----PARAMETERS-----%

prog = 5000;
E = 10;         % Epsilon parameter in the hamiltonian
temp = [0, 0];  % Temperature for simulation
flag_p = false; % Flag for doing plots
iter_max = 20;  % Maximum number of trials to find a valid config

%-----LOAD PARAMETERS-----%
if exist('params', 'var')
    parameters = who;
    for p = 1:length( parameters )
        if isfield(params, parameters{p})
            eval([parameters{p},'=','params.', parameters{p},';']); 
        end
    end
end

%-----PARAMETERS-----% END






%% %-----INITIALIZATION-----%

% Connectome size
N = size(C,1);
n = N;
D = sum(C(:))/ (N^2);

% Time steps
T = round(5*D*(N^4));

% Move probabilities ( jump (1-D), swap(D-1/N), flip(1/N) ) (cumsum)
R = [ 1-D, max(1-D,1-1/N) ];

% Hamiltonian tracking
stepsize = ceil( T / prog );
H = zeros(1, ceil(T/stepsize) );

% Function to evaluate transition probability
if all( temp == 0 )
    tran = @(a,b) 1*(b<0);
else
    tran = @(a,b) exp( -1 * b / (temp(1) + ((a-2)/(T-2))*(temp(2)-temp(1)) ) );
end

% Parent connectome barcode-list generation
[K(:,1), K(:,2), ~] = find( C );
B = size(K,1);

% (Barcode-pair,synapse) list, initialized randomly
S = randi(N, B, 2);

% Barcode matrix, (barcode,cell) and populate it from S
M = zeros(n,N);
for b = 1:B
    for p = 1:2
        M(K(b,p),S(b,p)) = M(K(b,p),S(b,p)) + 1;
    end
end

% Number of barcodes per cell
D = sum(M);

% Energy function
L = ...
    sum( -(1+E)*sum(M.^2) + E*(sum(M).^2) ) + ...   % Current E contrib.
    sum( sum(C+C').^2 );                % Normalize GS
H(1) = L;

%-----INITIALIZATION-----% END






%% %-----SIMULATION-----%

for t = 1:T
    r1 = rand;
    r2 = rand;
    
    %-----RECORD-----% BEGIN
    if stepsize == 1
        H( t ) = L;
    elseif mod(t, stepsize) == 1
        H( ceil(t/stepsize) ) = L;
    end
    %-----RECORD-----% END
    
    %-----ACTION-----% BEGIN
    if r1 < R(1)
        %-----JUMP-----% BEGIN
        b = ceil( B * rand );   % Select a barcode
        p = ceil( 2 * rand );   % Select an end
        c1 = S(b,p);            % Changing cell is S(b,p)
        c2 = c1;                % Changed cell is a random cell not c1
        while c2 == c1
            c2 = ceil( N * rand );
        end
        m = K(b,p);             % The barcode that moves is K(b,p)
        % Change in energy
        dE = 2 * ( E*(1+D(c2)-D(c1)) - (1+E)*(1+M(m,c2)-M(m,c1)) );
        
        % Accept clause
        if r2 < tran(t,dE)
            D(c1) = D(c1) - 1;
            D(c2) = D(c2) + 1;
            M(m,c1) = M(m,c1) - 1;
            M(m,c2) = M(m,c2) + 1;
            S(b,p) = c2;
            L = L + dE;
        end
        %-----JUMP-----% END
    elseif r1 < R(2)
        %-----SWAP-----% BEGIN
        p = 0;                  % Use p to count the recurrence
        b = [];                 % Use b as an index storage
        while ( length(b) < 2 ) && ( p <= iter_max )
            p = p + 1;
            c = ceil( N * rand );
            b = find( any(S==c,2) );
        end
        if p > iter_max         % If the last loop iterated a lot, stop
            continue
        end
        % Choose randomly those two barcode-pair ID's from the cell
        b1 = ceil( length(b) * rand );
        b2 = b1;
        while b2 == b1
            b2 = ceil( length(b) * rand );
        end
        b1 = b(b1);
        b2 = b(b2);
        % Get the pre-pos ends that the barcode-pair point to
        if S(b1,1) ~= c
            p1 = 1;
        elseif S(b1,2) ~= c
            p1 = 2;
        else % Pick random end if a self synapse of the cell
            p1 = ceil( 2 * rand );
        end
        if S(b2,1) ~= c
            p2 = 1;
        elseif S(b2,2) ~= c
            p2 = 2;
        else
            p2 = ceil( 2 * rand );
        end
        c1 = S(b1,p1);                  % Assign cells
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
        if r2 < tran(t,dE)
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
        b = ceil( B * rand );           % Randomly choose a barcode-pair
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
        if r2 < tran(t,dE)
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
        T = t;
        break
    end
    %-----TRUNCATE-----% END
end

%-----SIMULATION-----% END






%% %-----RECONSTRUCTION-----%

% Assign each cell a number that indicates their barcode ID
% Identify cells by most prominant barcode
conid = isOBOCtest( M );
% check for repeating barcode-id's
if conid == false
    success = 0;    % Fail if not OBOC
else
    C2 = sparse( conid(S(:,1)), conid(S(:,2)), 1, max(N,n), max(N,n));
    if all( C(:)==C2(:) )
        success = T;
    else
        success = 0;
    end
end
%-----RECONSTRUCTION-----% END






%% %-----PLOTTING-----%
if flag_p && (~isempty(conn{2}))
    f = figure(2);
    f.WindowStyle = 'docked';
    I = zeros(N,N,3);
    I(:,:,1) = C;
    I(:,:,2) = C2;
    subplot(1,2,1);
    imagesc(I); axis image;
    title(['Reconstruction (in ', num2str(T), ' moves)']);
    xlabel('Pre-synaptic cell ID');
    ylabel('Post-synaptic cell ID');
    subplot(1,2,2);
    semilogy( H( 1:find(H==0,1) ) );
    title('Hamiltonian evolution');
    if stepsize == 1
        xlabel('Steps');
    else
        xlabel(['Steps (x', num2str(stepsize), ')']);
    end
    ylabel('Energy above GS')
end

%-----PLOTTING-----% END

end






%% %-----AUXILLARY-----%
function [ map ] = isOBOCtest( marker )
% ISOBOC Returns cell-ID map if barcodes are OBOC, and false otherwise.
%   mapping is such that answer(cell_number) will give the cell as
%   identified by barcode.
%   Requires the marker matrix
% This makes it so that both cases where;
%   -Parent connectome is smaller than the tabula rasa one
% 	-Parent connectome is the same size as the tabula rasa one, but has 
%   unconnected cells
%   -Parent connectome is larger than tabula-rasa, but has empty cells to 
%   accomodate an OBOC solution on tabula-rasa network
%       -> If the parent is larger, but connected is as much as child, no
%       emptiness will occur
%       -> If the parent is same size, but connected is less than child;
%       some barcodes will be empty and they will show up in the difference
%       betwee {1,2,...,child}.
%       -> If the parent is larger but connected is less than child;
%           -> Worst case is all connected ID's are less than child_size.
%           Then unconnected will appear within the range [1, child]
%           sufficiently that random assigning will handle it
%           -> If disconnected barcodes have larger ID# than child_size,
%           they will increase the number of available free assigned ID's.
%           The difference [1,...,child_size]\[id'd_barcodes] will be
%           larger than empty cells in child
% Marker matrix of form (barcode ID, cell#)

% Initialize
map = false;

% Colsum should be equal to colmax
if any( max(marker) ~= sum(marker) )
    % Ensures that all barcodes in a cell is a single barcode type
    return
% Rowsum should be equal to rowmax
elseif any( max(marker,[],2) ~= sum(marker,2) )
    % Ensures that a barcode type is only present in a single cell
    return
end

% Not interested in how many barcodes now, but which barcode ID is present
[~,map] = max( marker );
% Remove the false 1 ID's of empty cells by O'ing their identification ID
map( sum(marker) == 0 ) = 0;
% Check which ones are missing ID's
miss = find( map==0 );
% Usable ID numbers for identification 
uuid = setdiff( 1:size(marker,2), map(map~=0) );
% Identify empty cells in order
map( miss ) = uuid( 1:length(miss) );

end
%-----AUXILLARY-----% END















































