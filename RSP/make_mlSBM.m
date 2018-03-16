function [A, P] = make_mlSBM(n,k,m,scaling_type,cin,lambda)
% n is the number of nodes
% k is the number of communities
% scaling_type is either constant or logarithmic

    cout = (1-lambda)*cin;
    
    for layers = 1:m
        seed = randi(1000);
        [G,P] = sbm_gen(n,k,cin,cout,seed);

        % Check for isolated nodes and fix
        isol = find(sum(G, 2) == 0);
        if ~isempty(isol)
            for idx = 1:numel(isol)
               neighbors = find(P == P(isol(idx)));
               edge = datasample(neighbors, 1);
               G(isol(idx),edge) = 1;
               G(edge, isol(idx)) = 1;
            end
        end

        %G(G == 2) = 1;
        A(:,:,layers) = full(G);
    end


end