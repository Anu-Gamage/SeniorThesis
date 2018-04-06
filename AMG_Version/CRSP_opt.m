function [acc_arr, nmi_arr, final_labels] = CRSP_opt(A, C, labels, n, k, m, b)
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
% Anuththari Gamage, 3/22/2018
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
    %P_joint = combine_P(P_ref);           % Combines probability matrix
    P_joint = combine_P_opt(P_ref);      % Combines probability matrix -- Xiaozhe
    P_joint = stochastize(P_joint);
    W = P_joint.*exp(-b*C_joint);              % Combined weights
    
    % no need to check in this way -- Xiaozhe
    %specRadius = abs(eigs(W,1, 'LM'));                    % Convergence check
    specRadius= max(abs(sum(W,2)));  % use the factor that W is nonnegative and apply PF theorem
    if specRadius >= 1
        error('Will not converge')
    end
    
    IW = speye(n) - W;
    Z = inv(IW);
    S = (Z*(C_joint.*W)*Z)./Z;
    C_bar = S - ones(n,1)*diag(S)';
    dRSP = (C_bar + C_bar')./2;   
    infFlag = 1e12;
    dRSP(isnan(dRSP)) = infFlag;                 % Flagging inf
    
    % Spectral Clustering   
%     aff = 1./(eye(n) + dRSP) - eye(n);      % Affinity Matrix
%     D = diag(sum(aff,2)) ;
%     L = (D^(-1/2))*aff*(D^(-1/2));          % Normalized Laplacian
%     [V,E] = eig(L);
%     [~,I] = sort(diag(E),'descend');
%     V = V(:, I(1:k)');      % the first one useful?
%     %V = V./sqrt(sum(V.^2,2)); % does not work on old version MATLAB
%     V=spdiags(1./sqrt(sum(V.^2,2)),0,n,n)*V;
    
    %-------------------------------
    % Sparse Spectral Clustering -- Xiaozhe
    %-------------------------------
    % find k-nearst-neighbors
    kn = 10;  % change this to change the knn graph
    [col_idx, val] = knnsearch((1:n)', (1:n)', 'K', kn+1, 'Distance', @(I,J)dist_CRSP(I,J,dRSP));
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
    amgParam.coarsest_size = max(10*k, 50);
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

function new_C = combine_C(C)
    m = size(C,2); 
    new_C = C{1};
    nz_C = C{1}~=0;         % Tracks count of non-zero costs
    for layers = 2:m
        new_C = new_C + C{layers};
        nz_C = nz_C + (C{layers}~=0);
    end
    %new_C = new_C./(nz_C + (nz_C==0));  % this is O(n^2) complexity
    % because nz_C + (nz_C==0) is basically a dense matrix -- Xiaozhe
    %-----------------
    %-- This approach only uses nonzeros, which is still O(n) -- Xiaozhe
    [row, col, val] = find(new_C); 
    [~, ~, val_nz] = find(nz_C);
    new_C = sparse(row, col, val./val_nz);
    %-----------------
end

function new_P = combine_P(P) %this subroutine requires O(n^2) operations -- Xiaozhe
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

function new_P = matrix_mult(P)  %  this subroutine requires O(n^2) operations -- Xiaozhe
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


%---------------------------
% combine P -- Xiaozhe
%---------------------------
function new_P = combine_P_opt(P) % if sparsity A + B is O(n), this is O(n) -- Xiaozhe
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

function P = stochastize(A)
    n = size(A,1);
   while sum(round(sum(A,2))' == ones(1,n)) ~= n
        %A = A./sum(A,1);
        %A = A./sum(A,2);
        %for i = 1:n                         % To avoid ./ issues (remove later)
        %     A(:,i) = A(:,i)/sum(A(:,i));
        %end
        A =  A*spdiags(1./sum(A,1)',0,n,n);  % avoid loop -- Xiaozhe
        %for i = 1:n
        %    A(i,:) = A(i,:)/sum(A(i,:));
        %end
        A = spdiags(1./sum(A,2),0,n,n)*A;  % avoid loop - Xiaozhe
   end
    P = A;
   %disp(sprintf('RowSums = %d',sum(round(sum(A,2))' == ones(1,n))))
end

% distance function for knnsearch -- Xiaozhe
function [ D ] = dist_CRSP(I, J, dRSP ) 

    D = dRSP(J,I);

end
