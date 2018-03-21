% Common - Randomized Shortest Paths
% Anuththari Gamage
% 2/17/2018

function [acc, nmiVal, final_labels] = CRSP(A, C, labels, n, k, m, b)
    
    % C-RSP
    Pref = zeros(n,n,m);
    for i = 1:m
        node_degrees = sum(A(:,:,i),2);
        x  = find(node_degrees);
        node_degrees(x) = 1./node_degrees(x);
        invD = diag(node_degrees);        % Inverse of Degree matrix
        Pref(:,:,i) = invD*A(:,:,i);      % Reference Transition Probability 
    end
       
    % Construct common W
    C_sum = combineP(C);
    %C_sum(C_sum==0) = inf;
    %C_sum = (sum(C,3));                     % Common cost matrix - normalize later
    
    combinedP = combineP(Pref);
    combinedP = stochastize(combinedP);
    W = combinedP.*exp(-b*C_sum);             % Combined weights
    
    specRadius = max(abs(eig(W)));
    if specRadius >= 1
        error('Will not converge')
    end
    
    Z = inv(eye(n) - W);
    S = (Z*(C_sum.*W)*Z)./Z;
    C_bar = S - ones(n,1)*diag(S)';
    dRSP = (C_bar + C_bar')./2;
    dRSP(isnan(dRSP)) = 1e6;                 % Flagging inf
    %disp(dRSP);
    %figure;imagesc(dRSP);colorbar 
    
    % Spectral Clustering
    aff = 1./(eye(n) + dRSP);               % Affinity Matrix
    %figure;imagesc(aff);colorbar
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
%     disp(labels)
%     disp(final_labels)

    % Accuracy
    acc = 100*sum(labels == final_labels)/n;
   % fprintf('Accuracy: %.2f\n', acc)
    nmiVal = nmi(labels, final_labels);
   % fprintf('NMI: %.2f\n', nmiVal)
end

function newP = combineP(P)
   n = size(P,1);
   %m = size(P,3);
   nz = sum(P~=0, 3); 
   Pnz = P + (P==0);
   
   newP = ones(n);
   for i = 1:size(P,3)
    newP = newP.*Pnz(:,:,i);
   end
   
   [r,c] = find(nz == 0);
   for i = 1:numel(r)
      newP(r(i),c(i)) = 0;
   end
      
   nz = nz + (nz == 0);         % To eliminate 0th roots
   newP = arrayfun(@(P,nz) nthroot(P,nz), newP, nz);
   

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