function [ connectome, I, U ] = cf_flickerSynapses( N, S, U )
% CF_FLICKERSYNAPSES: Generate active synapses list

% Number of active synapses requested
U = round( N * N * U );

% Activate symmetric map for the new random synapses
connectome = triu(rand(N));
sto = sort( connectome(connectome(:)~=0) );

connectome = ( connectome + connectome' + diag(rand(N,1)) );
connectome = connectome <= sto(U);

% Reactivate occupied synapses
for b = 1:size(S,1)
    connectome( S(b,1), S(b,2) ) = true;
    connectome( S(b,2), S(b,1) ) = true;
end

% Generate I and U
I = zeros(N);
U = zeros(N,1);
for b = 1:N
    sto = find( connectome(b,:);
    U(b) = length(sto);
    I(b, 1:U(b)) = sto;
end

end