% Randomized Shortest Paths
% Anuththari Gamage
% 2/17/2018

function [acc, nmiVal, final_labels] = RSP(A, labels, n, k, b)
    
    % RSP
    node_degrees = sum(A,2);
    x  = find(node_degrees);
    node_degrees(x) = 1./node_degrees(x);
    invD = diag(node_degrees);        % Inverse of Degree matrix
    Pref = invD*A;                    % Reference Transition Probability 

    C = A;                            % Cost matrix (defined to have equal costs here)
    W = Pref.*expm(-b*C);
    specRadius = max(abs(eig(W)));
    if specRadius >= 1
        error('Will not converge')
    end
    
   Z = inv(eye(n) - W);
   S = (Z*(C.*W)*Z)./Z;
   C_bar = S - ones(n,1)*diag(S)';
   dRSP = (C_bar + C_bar')./2;
   dRSP(isnan(dRSP)) = -1;                 % Flagging inf
   %disp(dRSP);
   %figure;imagesc(dRSP);colorbar 


    % Spectral Clustering
    aff = 1./(eye(n) + dRSP);               % Affinity Matrix
  %  figure;imagesc(aff);colorbar 
    L = diag(n*ones(n,1)) - aff;
    [V,~] = eig(L);
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
    nmiVal = nmi(labels, final_labels);

end