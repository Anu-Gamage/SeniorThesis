% Main test file for C-RSP vs Everything else
% Anuththari Gamage
% 3/24/2018
clear;clc;close all

n = 500;                        % no. of nodes 
k = 2;                          % no. of clusters
m_array = [2];                % no. of layers
b = 0.02;                       % Tuning parameter for CRSP
c = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0,15.0, 20.0];               % Varying node degree
lambda = 0.9;                   % For sbm-gen
lambda_scml = 0.5;              % regularization parameter for SC-ML
lambda_coreg = 0.01;            % Co-regularization parameter for CPSC/CCSC
num_iter = 10;                  % no. of iterations for CPSC
do_plot = 0;                    % To plot data matrices
do_result_plot = 1;             % To plot results
num_runs = 2;                  % Number of runs > 1

% Index of algorithms used for comparison: 
% 1=C-RSP, 2=SC-ML, 3=C-FE, 4=CPSC, 5=CCSC, 6=MultiNMF, 7= CRSP_opt, 8=CFE_opt
num_compar = [1,6];                 
alg_names = [string('CRSP'), string('SCML'), string('CFE'), string('CPSC'),string('CCSC'), string('MultiNMF'), string('CRSP-opt'), string('CFE-opt')];
ccr_array = zeros(num_runs, numel(c), numel(m_array), numel(num_compar));
nmi_array = zeros(num_runs, numel(c), numel(m_array), numel(num_compar));


for runs = 1:num_runs
    % Generate adjacency tensor of test data
    data = cell(1,numel(c));   
    for i = 1:numel(c)
        [data{i}, labels] = mlsbm_gen(n,k,max(m_array), c(i), lambda);
    end
    for alg = 1:numel(num_compar)
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
                switch num_compar(alg)
                    case 1
                    [acc(degree), nmi(degree), final_labels] = CRSP(A,A,labels', n, k,m, b); % Cost matrix = A
                    disp('C-RSP')
                    fprintf('CCR: %.2f\n', acc(degree))
                    fprintf('NMI: %.2f\n\n', nmi(degree))
                    case 2
                    addpath([pwd '/SCML/']) 
                    disp('SC-ML')
                    [acc(degree), nmi(degree), final_labels] = SCML(A,k,lambda_scml, labels);
                    fprintf('CCR: %.2f\n', acc(degree))
                    fprintf('NMI: %.2f\n\n', nmi(degree))
                    case 3
                    [acc(degree), nmi(degree), final_labels] = CFE(A,A,labels', n, k,m, b); % Cost matrix = A
                    disp('C-FE')
                    fprintf('CCR: %.2f\n', acc(degree))
                    fprintf('NMI: %.2f\n\n', nmi(degree))
                    case 4
                    sigma = zeros(1,m);
                    addpath([pwd '/coregularizedSC/'])  
                    [acc(degree), nmi(degree),final_labels] = spectral_pairwise_multiview(A,m,k,sigma,lambda_coreg, labels', num_iter);          
                    disp('CPSC')
                    fprintf('CCR: %.2f\n', acc(degree))
                    fprintf('NMI: %.2f\n\n', nmi(degree))
                    case 5
                    sigma = zeros(1,m);
                    addpath([pwd '/coregularizedSC/']) 
                    [acc(degree), nmi(degree),final_labels] = spectral_centroid_multiview(A,m,k,sigma,repmat(lambda_coreg,1,m), labels, num_iter);          
                    disp('CCSC')
                    fprintf('CCR: %.2f\n', acc(degree))
                    fprintf('NMI: %.2f\n\n', nmi(degree))
                    case 6
                    addpath([pwd '/MultiNMF/']) 
                    disp('MultiNMF')
                    options = [];               % Default values below
                    options.maxIter = 75;      % 100
                    options.error = 1e-6;       % 1e-6
                    options.nRepeat = 30;       % 30
                    options.minIter = 50;       % 50
                    options.meanFitRatio = 0.1; % 0.1
                    options.rounds = 20;        % 20
                    options.WeightMode='Binary';
                    options.varWeight = 0;      % 1
                    options.kmeans = 1;         % 1
                    options.Gaplpha= 0;        	% 1     graph regularization parameter                      
                    options.alpha=0.;          % 0.1
                    options.delta = 0.1;        % 0.1
                    options.beta = 0;          % 10
                    options.gamma = 2;          % 2
                    options.K = k;
                    options.alphas = ones(1,m);         % Equal weights to each layer
                    
                    [acc(degree), nmi(degree)] = GMultiNMF(A, k, A, labels, options);                    
                    fprintf('CCR: %.2f\n', acc(degree))
                    fprintf('NMI: %.2f\n\n', nmi(degree))
                    
                    case 7
                    addpath([pwd '/AMG_version/']) 
                    addpath([pwd '/AMG_version/Matlab-AMG'])
                    [acc(degree), nmi(degree), final_labels] = CRSP_opt(A,A,labels', n, k,m, b); % Cost matrix = A
                    disp('C-RSP_opt')
                    fprintf('CCR: %.2f\n', acc(degree))
                    fprintf('NMI: %.2f\n\n', nmi(degree))
                    
                    case 8
                    addpath([pwd '/AMG_version/']) 
                    addpath([pwd '/AMG_version/Matlab-AMG'])
                    [acc(degree), nmi(degree), final_labels] = CFE_opt_approx(A,A,labels', n, k,m, b); % Cost matrix = A
                    disp('C-FE_opt')
                    fprintf('CCR: %.2f\n', acc(degree))
                    fprintf('NMI: %.2f\n\n', nmi(degree))
                end
            end
            ccr_array(runs,:, layers, alg) = acc;
            nmi_array(runs,:, layers, alg) = nmi;
        end
    end
end
% save('ccr.mat', 'ccr_array')
% save('nmi.mat', 'nmi_array')

if do_result_plot   
     x =  [ 0    0.4470    0.7410;    0.8500    0.3250    0.0980;    0.9290    0.6940    0.1250; 0.4940    0.1840    0.5560; 0.4660    0.6740    0.1880;  0.3010    0.7450    0.9330;  0.6350    0.0780    0.1840]; 
     for i = 1:numel(m_array) 
         figure
         for alg = 1:numel(num_compar)
            avg_ccr = mean(ccr_array(:,:,i, alg));
            avg_nmi = mean(nmi_array(:,:,i, alg));
            std_ccr = std(ccr_array(:,:,i, alg));
            std_nmi = std(nmi_array(:,:,i, alg));
            yyaxis left; errorbar(c, avg_ccr, std_ccr, '-', 'color', x(alg,:),'DisplayName', sprintf('%s-CCR', alg_names(num_compar(alg)))); 
            ylabel('CCR'); ylim([ 50,100])
%             hold on;yyaxis right; errorbar(c, avg_nmi, std_nmi, '--', 'color', x(alg,:),'DisplayName',sprintf('%s-NMI', alg_names(num_compar(alg)))); 
%             ylabel('NMI'); ylim([0,1])
            hold on
         end
        
        title_string = alg_names(num_compar(1));
        for idx = 2:numel(num_compar)
            title_string = strcat(title_string, string(' vs. '), alg_names(num_compar(idx)));
        end
        title(sprintf('%s: Nodes = %d, Clusters = %d, Layers = %d, Runs = %d',title_string, n, k,m_array(i), num_runs))
        xlabel('c'); 
        legend('Location','SouthEast')  
        
        title_string = alg_names(num_compar(1));
        for idx = 2:numel(num_compar)
            title_string = strcat(title_string, string('_'), alg_names(num_compar(idx)));
        end
       % saveas(gcf, [pwd '/figs/' sprintf('%s_n%d_k%d_m%d_r%d.png', title_string,n,k,i,num_runs)])
     end
end