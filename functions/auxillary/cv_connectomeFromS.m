function [ connect ] = cv_connectomeFromS ( synapseList, N )
% CONNFROMS generates connectivity (logical) from synapse matrix
    synapseList = unique(synapseList, 'rows');
    connect = full(sparse(synapseList(:,1),synapseList(:,2),true,N,N));
end