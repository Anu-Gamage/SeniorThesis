% Main test file for C-RSP
% Anuththari Gamage
% 3/16/2018
clear;clc;close all

n = 500;                        % no. of nodes 
k = 2;                          % no. of clusters
m_array = [1,2,3];              % no. of layers
b = 0.02;                       % Tuning parameter for CRSP
c = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0,15.0, 20.0];               % Varying node degree
lambda = 0.9;                   % For sbm-gen
lambda_scml = 0.5;              % regularization parameter for SC-ML
do_plot = 0;                     % To plot data matrices
do_result_plot = 1;              % To plot results
num_runs = 15;

num_compar = 2;                 % no. of algorithms used for comparison
ccr_array = zeros(num_runs, numel(c), numel(m_array), num_compar);
nmi_array = zeros(num_runs, numel(c), numel(m_array), num_compar);


for runs = 1:num_runs
    % Generate adjacency tensor of test data
    data = cell(1,numel(c));   
    for i = 1:numel(c)
        [data{i}, labels] = mlsbm_gen(n,k,max(m_array), c(i), lambda);
    end
    for alg = 1:num_compar
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
                % Run algorithms
                switch alg
                    case 1
                    [acc(degree), nmi(degree), final_labels] = CRSP(A,A,labels', n, k,m, b); % Cost matrix = A
                    disp('C-RSP')
                    fprintf('CCR: %.2f\n', acc(degree))
                    fprintf('NMI: %.2f\n\n', nmi(degree))
                    case 2
                    disp('SC-ML')
                    [acc(degree), nmi(degree), final_labels] = SCML(A,k,lambda_scml, labels);
                    fprintf('CCR: %.2f\n', acc(degree))
                    fprintf('NMI: %.2f\n\n', nmi(degree))
                end
            end
            ccr_array(runs,:, layers, alg) = acc;
            nmi_array(runs,:, layers, alg) = nmi;
        end
    end
end
% save('crsp_ccr.mat', 'ccr_array')
% save('crsp_nmi.mat', 'nmi_array')

if do_result_plot   
     for i = 1:numel(m_array) 
         figure
         for alg = 1:num_compar
            avg_ccr = mean(ccr_array(:,:,i, alg));
            avg_nmi = mean(nmi_array(:,:,i, alg));
            std_ccr = std(ccr_array(:,:,i, alg));
            std_nmi = std(nmi_array(:,:,i, alg));
            yyaxis left; errorbar(c, avg_ccr, std_ccr); ylabel('CCR'); ylim([ 50,100])
            hold on;yyaxis right; errorbar(c, avg_nmi, std_nmi); ylabel('NMI'); ylim([0,1])
            hold on
         end
        title(sprintf('CRSP vs SC-ML: Nodes = %d, Clusters = %d, Layers = %d', n, k,m_array(i)))
        xlabel('c'); 
        legend('CRSP-CCR', 'CRSP-NMI','SCML-CCR', 'SCML-NMI','Location','SouthEast')
     end
end