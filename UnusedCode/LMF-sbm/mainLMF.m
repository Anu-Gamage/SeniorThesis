% Main test file for LMF
% Anuththari Gamage
% 3/18/2018
clear;clc;close all

n = 100;                        % no. of nodes 
k = 2;                          % no. of clusters
m_array = [1];                  % no. of layers
alpha = 1e-5;                   % Correction parameter for LMF
c = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0,15.0, 20.0];               % Varying node degree
lambda = 0.9;
do_plot = 0;                     % To plot data matrices
do_result_plot = 1;              % To plot results
num_runs = 5;
dispError = 0;
error_threshold = 1e-5;          % Error threshold for optimization
ccr_array = zeros(num_runs, numel(c), numel(m_array));
nmi_array = zeros(num_runs, numel(c), numel(m_array));


for runs = 1:num_runs
    % Generate adjacency tensor of test data
    data = cell(1,numel(c));   
    for i = 1:numel(c)
        [data{i}, labels] = mlsbm_gen(n,k,max(m_array), c(i), lambda);
    end
    for layers = 1:numel(m_array)        % Varying no. of layers
        m = m_array(layers);
        acc = zeros(1, numel(c));
        nmi = zeros(1, numel(c));
        for degree = 1:numel(c)         % Varying node degree
            A = cell(1,layers);
            for i = 1:m                 % Select relevant tensor
                A{i} = data{degree}{i};
            end
            fprintf('Variable %d processing:\n', degree)
            if do_plot
               figure; 
               for i = 1:m
                  subplot(1,m,i);spy(A{i}); title(sprintf('Layer %d', i))
               end
            end
            % Run LMF
            for i = 1:m
                G(:,:,i) = full(A{i});
            end
            A = G;
            [acc(degree), nmi(degree), final_labels] = lmf(A,labels',k,alpha,dispError, error_threshold); 
            fprintf('CCR: %.2f\n', acc(degree))
            fprintf('NMI: %.2f\n\n', nmi(degree))
        end
        ccr_array(runs,:, layers) = acc;
        nmi_array(runs,:, layers) = nmi;
    end
end
save('lmf_ccr.mat', 'ccr_array')
save('lmf_nmi.mat', 'nmi_array')

if do_result_plot
    for i = 1:numel(m_array)
            avg_ccr = mean(ccr_array(:,:,i));
            avg_nmi = mean(nmi_array(:,:,i));
            std_ccr = std(ccr_array(:,:,i));
            std_nmi = std(nmi_array(:,:,i));
            figure;yyaxis left; errorbar(c, avg_ccr, std_ccr); ylabel('CCR'); ylim([ 50,100])
            hold on;yyaxis right; errorbar(c, avg_nmi, std_nmi); ylabel('NMI'); ylim([0,1])
            title(sprintf('LMF: Nodes = %d, Clusters = %d, Layers = %d', n, k,m_array(i)))
            xlabel('c'); 
            legend('CCR', 'NMI','Location','SouthEast')
    end
end