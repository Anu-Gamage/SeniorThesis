function [acc_val, nmi_val, final_labels] = SCML(A, k, lambda_scml, labels)
% Modified function file for SC-ML
% Inputs:
% A - multilayer graph tensor (full, not sparse)
% labels - true cluster assignments of nodes
% k - no. of clusters
% Outputs:
% acc_arr - accuracy 
% nmi_arr - normalized mutual information
% final_labels - estimated labels from SC-ML
% Anuththari Gamage, 3/19/2018

        % Convert to dense matrix
        n = size(A{1},1);
        m = size(A,2);
        G = zeros(n,n,m);
        for i = 1:m
            G(:,:,i) = full(A{i});
        end
        A = G;    
        
        % Run SC-ML
        final_labels = sc_ml(A,k,lambda_scml); 

        % Re-order labels
        final_labels = orderLabels(final_labels, labels);

        % Accuracy
        acc_val = 100*sum(labels == final_labels)/n;
        nmi_val = nmi(labels', final_labels');       
end

function final_labels = orderLabels(est_labels, labels)
    
    %Linear Sum Assignment
    final_labels = zeros(numel(est_labels),1);

    C = confusionmat(labels, est_labels);
    new_labels = munkres(-C);
    if ~isequal(new_labels, 1:length(new_labels))
       for i = 1:length(new_labels)
           final_labels(est_labels == new_labels(i)) = i;
       end
    else
        final_labels = est_labels;
    end
end