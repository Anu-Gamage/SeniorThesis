function [A, labels] = mlsbm_gen(n,k,m,cin,lambda)
% Generates multi-layers SBM graphs
% Inputs:
% n - number of nodes
% k - number of communities
% m - number of layers
% cin - intra-cluster node degree
% lambda - intra-cluster node degree
% Outputs: 
% A - cell array of SBM graphs
% labels - cluster assignment for each node
% Anuththari Gamage, 3/16/2018 

    A = cell(1,m);
    cout = (1-lambda)*cin;
    
    for layers = 1:m
        seed = randi(1000);
        %[G,labels] = sbm_gen(n,k,cin,cout,seed);
        perc = [0.3,0.7];
        [G,labels] = unbalanced_sbm_gen(n,2,cin,cout,seed, perc);

        % Check for isolated nodes and fix
        isol = find(sum(G, 2) == 0);
        if ~isempty(isol)
            for idx = 1:numel(isol)
               neighbors = find(labels == labels(isol(idx)));
               edge = datasample(neighbors, 1);
               G(isol(idx),edge) = 1;
               G(edge, isol(idx)) = 1;
            end
        end

        A{layers} = G;
    end
end
