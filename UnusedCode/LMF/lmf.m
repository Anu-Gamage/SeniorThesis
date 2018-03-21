% Linked Matrix Factorization
% Anuththari Gamage
% 1/14/2018


function [acc, nmiVal, final_labels] = lmf(data, labels, num_clusters, alphaVal, dispError)

    global n k m A alpha P D
    A = data;
    n = size(A,1);
    m = size(A, 3);
    k = num_clusters;
    alpha = alphaVal;

    % Optimize P and D
    D = rand(k,k,m);
    P = rand(n,k);
    fvalP_old = 0;
    fvalD_old = 0;
    fvalP = 10;
    fvalD = 10;
    %numIter = 20;
    iter = 0;

    while((fvalP ~= fvalP_old && fvalD ~= fvalD_old) && iter < 1500)
        fvalP_old = fvalP;
        fvalD_old = fvalD;

        iter = iter + 1;
    %     disp(iter)

        % Generate P fixing D
        P = reshape(P, [1,size(P,1)*size(P,2)]);
        options = struct('Display', 'off');
        [P,fvalP,exitflag,output] = lbfgs(@objFuncP, P, options);
        P = reshape(P, n,k);
        
        if dispError
            fprintf('Error for P: %f\n', fvalP)
        end

        % Generate D fixing P
        global j
        for j = 1:m
            D_m = reshape(D(:,:,j), [1, k*k]);
            options = struct('Display', 'off');
            [D_m,fvalD,exitflag,output] = lbfgs(@objFuncD, D_m, options);
            D(:,:,j) = reshape(D_m, [k,k]);
            
            if dispError
                fprintf('Error for D(%d): %f\n\n', j, fvalD)
            end
        end
    end

    % Clustering using k-Means and Linear Sum Assignment
    est_labels = kmeans(P, k)';
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

%%% FUNCTIONS %%%

function [error, grad] = objFuncP(P0)    % 1st arg - fval, 2nd arg - gradient vector same size as x
    global A m n k D alpha
    error = 0;
    P = reshape(P0, [n,k]);
    
    for i = 1:m
        error = error + 0.5*((norm(A(:,:,i) - P*D(:,:,i)*P', 'fro'))^2) + (alpha/2)*(norm(D(:,:,i),'fro')^2 + norm(P, 'fro')^2);
    end
    
    grad = 0;
    for i = 1:m
       grad = grad + (A(:,:,i)-P*D(:,:,i)*P')*P*D(:,:,i); 
    end
    grad = -2*grad + alpha*P;
    grad = reshape(grad, [1, size(grad,1)*size(grad,2)]);
end


function [error, grad] = objFuncD(D0)    % 1st arg - fval, 2nd arg - gradient vector same size as x
    global A k P alpha j
    
    D_m = reshape(D0, [k,k]);  
    error = 0.5*((norm(A(:,:,j) - P*D_m*P', 'fro'))^2) + (alpha/2)*(norm(D_m,'fro')^2 + norm(P, 'fro')^2);
    grad = -P'*(A(:,:,j) - P*D_m*P')*P + alpha*D_m;
    grad = reshape(grad, [1, size(grad,1)*size(grad,2)]);
end

