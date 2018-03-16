% Common - Randomized Shortest Paths
% Inputs:
% A - multilayer graph tensor
% C - cost tensor
% labels - true cluster assignments of nodes
% n - no. of nodes
% k - no. of clusters
% m - no. of layers
% b - RSP tuning paramter
% Outputs:
% acc_arr - accuracy 
% nmi_arr - normalized mutual information
% final_labels - estimated labels from RSP
% Anuththari Gamage, 3/16/2018
function [acc_arr, nmi_arr, final_labels] = CRSP(A, C, labels, n, k, m, b)
    P_ref = cell(1,m);                  % Reference transition probability
    for i = 1:m
        node_degrees = sum(A{i},2);
        inv_D = sparse(1:n, 1:n, 1./node_degrees);  % Inverse of Degree matrix
        P_ref{i} = inv_D*A{i};      
    end
       
    % Construct common W
    C_joint = matrix_mult(C);                  % Combined cost matrix   
    P_joint = combine_P(P_ref);                % Combines probability matrix
    P_joint = stochastize(P_joint);
    W = P_joint.*exp(-b*C_joint);              % Combined weights
    
    specRadius = eigs(W,1);                    % Convergence check
    if specRadius >= 1
        error('Will not converge')
    end
    
    Z = inv(speye(n) - W);
    S = (Z*(C_joint.*W)*Z)./Z;
    C_bar = S - ones(n,1)*diag(S)';
    dRSP = (C_bar + C_bar')./2;
    dRSP(isnan(dRSP)) = 1e6;                 % Flagging inf
    
    % Spectral Clustering
    aff = 1./(eye(n) + dRSP);                % Affinity Matrix
    L = diag(n*ones(n,1)) - aff;
    [V,D] = eig(L);
    vec = V(:,1:k);

    % Clustering using k-Means and Linear Sum Assignment
    est_labels = kmeans(vec, k)';
    final_labels = zeros(1, n);

    C = confusionmat(labels, est_labels);
    new_labels = munkres(-C);
    if ~isequal(new_labels, 1:length(new_labels))
       for i = 1:length(new_labels)
           final_labels(est_labels == i) = new_labels(i);
       end
    else
        final_labels = est_labels;
    end

    % Accuracy
    acc_arr = 100*sum(labels == final_labels)/n;
    nmi_arr = nmi(labels, final_labels);
end

function new_P = combine_P(P)
   m = size(P,2);              
   new_P = matrix_mult(P);
   
   % Taking nth roots
   roots = (P{1}~=0);
   for i = 2:m
       roots = roots + (P{i}~=0);
   end
   roots = roots + (roots==0);          % To eliminate 0th roots
   new_P = nthroot(new_P, roots);
end

function new_P = matrix_mult(P)
   m = size(P,2); 
   % Multiplication retaining all edges
   no_edges = (P{1} == 0);              % Record entries with no edges
   for i = 2:m
       no_edges =  no_edges.*(P{i}==0);
   end    
   
   Pnz = cell(1,m);                 % Add 1 to avoid multiplication error
   for i = 1:m
       Pnz{i} = P{i} + (P{i}==0);    
   end
   new_P = Pnz{1};
   for i = 2:m
    new_P = new_P.*Pnz{i};
   end
   new_P = new_P.*(no_edges ~= 1);       % Zero out entries with no edges 
end


function P = stochastize(A)
    n = size(A,1);
   while sum(round(sum(A,2))' == ones(1,n)) ~= n
        A = A./sum(A,1);
        A = A./sum(A,2);
   end
    P = A;
   %disp(sprintf('RowSums = %d',sum(round(sum(A,2))' == ones(1,n))))
end