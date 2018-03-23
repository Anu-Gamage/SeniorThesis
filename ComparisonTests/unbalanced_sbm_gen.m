function [A,conf_true] = unbalanced_sbm_gen(N,q,cin,cout,seed, perc)
% returns adjacency matrix and community membership vector generated by the
% stochastic block model with N nodes, q communities, intra-community
% average connectivity cin, extra-community average connectivity cout.
% seed initializes the random number generator.
%
% Modified to generate unbalanced clusters
% perc - percentage of N nodes contained in each cluster e.g. [0.3 0.7]
% Anuththari Gamage 3/22/2018

rng(seed);
size = floor(perc*N);
conf_true = zeros(N,1);
start_id = [1 cumsum(size(1:end-1))+1];
end_id = cumsum(size);

for k = 1:q
    conf_true(start_id(k):end_id(k)) = k;
end
conf_true(conf_true==0) = 1;

A = sparse(N,N);
for i = 1:q
    current_block = sprand(size(i),size(i),cin/N);
    current_block = (current_block ~= 0);
    current_block = triu(current_block,1);
    current_block = current_block + current_block';
    % keyboard
    A(start_id(i):end_id(i),start_id(i):end_id(i)) = current_block;
    for j = i+1:q
        current_block = sprand(size(i),size(j),cout/N);
        current_block = (current_block ~= 0);
        A(start_id(i):end_id(i),start_id(j):end_id(j)) = current_block;
        A(start_id(j):end_id(j),start_id(i):end_id(i)) = current_block';
    end
end
end
