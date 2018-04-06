function [acc_arr, nmi_arr, final_labels] = CFE_opt_approx(A, C, labels, n, k, m, b)
% Common - Free Energy Distance
% Inputs:
% A - multilayer graph tensor
% C - cost tensor
% labels - true cluster assignments of nodes
% n - no. of nodes
% k - no. of clusters
% m - no. of layers
% b - FE tuning paramter
% Outputs:
% acc_arr - accuracy 
% nmi_arr - normalized mutual information
% final_labels - estimated labels from RSP
% Anuththari Gamage, 3/24/2018
% Optimized and added AMG: Xiaozhe Hu, 4/5/2018

    P_ref = cell(1,m);                  % Reference transition probability
    for i = 1:m
        node_degrees = sum(A{i},2);
        %inv_D = sparse(1:n, 1:n, 1./node_degrees);  % Inverse of Degree matrix
        inv_D = spdiags(1./node_degrees, 0, n, n); % use spdiags to save memory -- Xiaozhe
        P_ref{i} = inv_D*A{i};      
    end
       
    % Construct common W
    C_joint = combine_C(C);                  % Combined cost matrix   
    P_joint = combine_P(P_ref);                % Combines probability matrix
    W = P_joint.*exp(-b*C_joint);              % Combined weights
    
    specRadius= max(abs(sum(W,2)));  % use the factor that W is nonnegative and apply PF theorem
    if specRadius >= 1
        error('Will not converge')
    end
    
%     Z = inv(speye(n) - W);
%     Dh = diag(diag(Z));
%     Zh = Z*inv(Dh);
    Z = speye(n); 
    temp = W;
    for i=1:15
        Z = Z + temp;
        temp = temp*W;
    end
    invDh = spdiags(1./diag(Z),0,n,n);
    Zh = Z*invDh;
    
%     phi = (-1/b)*log(Zh);
%     dFE = (phi + phi')./2;   
%     infFlag = 1e12;
%     dFE(dFE == inf) = infFlag;                 % Flagging inf

    [row, col, val_Zh] = find(Zh);
    phi = sparse(row,col,(-1/b).*log(val_Zh),n,n);
    dFE = (phi + phi')./2;
    dFE = dFE - spdiags(diag(dFE),0,n,n);
    
    %-------------------------------
    % Sparse Spectral Clustering -- Xiaozhe
    %-------------------------------
%     % find k-nearst-neighbors
    kn = 10;  % change this to change the knn graph
    %[col_idx, val] = knnsearch((1:n)', (1:n)', 'K', kn+1, 'Distance', @(I,J)dist_CFE(I,J,dFE));
    [col_idx, val] = knnsearch((1:n)', (1:n)', 'K', kn+1, 'Distance', @(I,J)dist_CFE_sparse(I,J,dFE));
    col_idx = col_idx(:,2:end)';
    col_idx = col_idx(:);
    row_idx = repmat((1:n), kn, 1);
    row_idx = row_idx(:);
    val = val(:,2:end);
    val = val(:);
    aff = sparse(row_idx,col_idx,1./val,n,n);
    aff = aff+aff';

      % use MATLAB eigensolver
%     invRootD = spdiags(1./sqrt(sum(aff,2)),0,n,n);
%     NA = invRootD*aff*invRootD;          % Normalize affinity (NA)
%     NA = (NA+NA')/2;                             % make sure NA is symmetric
%     [V, ~] = eigs(NA, k+1, 'LA');
%     V = V(:,2:k+1);
    
    % cascadic algebraic multigrid solver
    D = spdiags(sum(aff,2),0,n,n);
    L = D - aff;  
    % AMG parameters
    [ amgParam ] = init_AMG_param;
    amgParam.print_level = -1; 
    amgParam.coarsest_size = max(10*k, 120);
    %amgParam.agg_type = 'MWM';
    amgParam.smooth_type = 'SGS';
    amgParam.n_postsmooth =2;
    amgParam.number_eigen = k;
    % AMG setup
    amgData = AMG_Setup(L, amgParam);
    % AMG solve
    V = cascade_eig(amgData, 1, amgParam);
    
    %V = V./sqrt(sum(V.^2,2)); % does not work on old version MATLAB
    V=spdiags(1./sqrt(sum(V.^2,2)),0,n,n)*V;
    %-------------------------------
    
    % Clustering using k-Means and Linear Sum Assignment
    est_labels = kmeans(V, k)';
    final_labels = zeros(1, n);

    C = confusionmat(labels, est_labels);
    new_labels = munkres(-C);
    if ~isequal(new_labels, 1:length(new_labels))
       for i = 1:length(new_labels)
           final_labels(est_labels == new_labels(i)) = i;
       end
    else
        final_labels = est_labels;
    end
    % Accuracy
    acc_arr = 100*sum(labels == final_labels)/n;
    nmi_arr = nmi(labels, final_labels);
end

%---------------------------
% combine C -- Xiaozhe
%---------------------------
function new_C = combine_C(C)
    m = size(C,2); 
    new_C = C{1};
    nz_C = C{1}~=0;         % Tracks count of non-zero costs
    for layers = 2:m
        new_C = new_C + C{layers};
        nz_C = nz_C + (C{layers}~=0);
    end
    [row, col, val] = find(new_C); 
    [~, ~, val_nz] = find(nz_C);
    new_C = sparse(row, col, val./val_nz);
end

%---------------------------
% combine P -- Xiaozhe
%---------------------------
function new_P = combine_P(P) % if sparsity A + B is O(n), this is O(n) -- Xiaozhe
   m = size(P,2);       
   
   % multiply p
   new_P = matrix_mult_opt(P);  
   
   % compute roots
   roots = spones(P{1});
   for i=2:m
       roots = roots + spones(P{i});
   end
   
   % take roots
   [row, col, val_p] = find(new_P);
   [~, ~, val_roots] = find(roots);
   val_p = nthroot(val_p,val_roots);
   
   % get new P
   new_P = sparse(row, col, val_p);
   
   % row stochastic
   new_P = new_P./(sum(new_P,2));  % Make row stochastic
end

%---------------------------
% multiply two matrix -- Xiaozhe
%---------------------------
function C = A_multi_B(A,B)  % if sparsity A + B is O(n), this is O(n) -- Xiaozhe

    C = (A+B - A.*spones(B) - B.*spones(A)) + (A.*B);

end
%---------------------------

%---------------------------
% multiply P -- Xiaozhe
%---------------------------
function new_P = matrix_mult_opt(P) % if sparsity sum_i P{i} is O(n), this is O(n) -- Xiaozhe
       m = size(P,2); 
       new_P = P{1};
       
       for i = 2:m
            new_P = A_multi_B(new_P, P{i});
       end
end
%---------------------------

function [ D ] = dist_CFE(I, J, dFE ) 

    D = dFE(J,I);

end

function [ D ] = dist_CFE_sparse(I, J, dFE ) 

    D = dFE(J,I);
    idx = (D==0);
    D(idx) = 1e12;
    D(I) = 0;

end
