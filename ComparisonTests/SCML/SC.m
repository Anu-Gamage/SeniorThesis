function [acc_arr, nmi_arr, final_labels] = SC(A,k,labels)
% Runs spectral clustering on sums of adjacency matrices
% Anuththari Gamage, 3/21/2018
    
    n = size(A{1},1);
    sumA = A{1};
    for i = 2:size(A,2);
        sumA = sumA + A{i};                % Sum of adjacency matrices
    end
    sumA = full(sumA);
    
    % Spectral Clustering   
    D = diag(sum(sumA,2)) ;
    L = (D^(-1/2))*sumA*(D^(-1/2));          % Normalized Laplacian
    %L = D - sumA;                             % Unnormalized Laplacian
    [V,E] = eig(L);
    [~,I] = sort(diag(E),'descend');            % k largest eigenvalues
    V = V(:, I(1:k)');
 %   V = V./sqrt(sum(V.^2,2));
      
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

