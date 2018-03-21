% Generate MLSBM graph
% Anuththari Gamage 
% 12/26/2017
%
% Generate an MLSBM type multilayer graph with the following paramters:
%   n = no. of nodes 
%   k = no. of clusters
%   m = no. of layers 
%   p = SNR parameter (2-3 for strong signal, 1 for weak)
%   density = node density (percentage)
function [A, true_labels] = mlsbm(n,k,m,p,density)

    % Assign clusters
    true_labels = zeros(1,n);
    for i = 1:k
       true_labels(i*floor(n/k) - floor(n/k) + 1:i*floor(n/k)) = i;
    end
    num  = sum(true_labels == 0);
    if(num > 0)
        true_labels(true_labels == 0) = randperm(k, num);
    end

    % Create Z matrix
    Z = zeros(n,k);
    for i = 1:k
        ind = find(true_labels == i);
        Z(ind, i) = 1;
    end

    % Create Adjacency tensor
    A = zeros(n,n,m);
    for i = 1:m
        % Create B matrix
        d = rand(1,k);
        x = triu(rand(k),1)./p;
        B = (diag(d) + x + x');
        A_m = Z*B*Z';

        % Create A layer
        probMat = ones(n);
        for v = 1:ceil(5*density*n*n)
            probMat(randi([1,n]),randi([1,n])) = rand();
        end
       % probMat = rand(n);
        A_m(A_m > probMat) = 1;
        A_m(A_m~=1) = 0;     
        
        A(:,:,i) = A_m;
%         figure; spy(A_m)
%         title(sprintf('A(%d)', i))
    end
    
    
    %disp('done')
end
