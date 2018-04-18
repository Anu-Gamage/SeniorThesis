% Main test file for C-RSP and others using external data (not SBM-generated)
% Anuththari Gamage
% 3/31/2018
clear;clc;close all

% Algorithm to run: 1=CRSP 2=CFE 3=SCML 4=CPSC 5=CCSC 6=MultiNMF
algs = [1,2,3];
alg_names = [string('CRSP'), string('CFE'), string('SCML'), string('CPSC'),string('CCSC'), string('MultiNMF')];

b = 1e-2;                          % RSP parameter: b = 10 for 3Sources, b=0.1 for multiview twitter
lambda_coreg = 0.01;             % Co-regularization parameter for CPSC/CCSC
num_iter = 10;                   % no. of iterations for CPSC
do_result_plot = 1;              % To plot results


% Data parameters
% addpath([pwd '/Datasets/MultiviewTwitter/'])  
% load('politicsuk_adjacency.mat')
% load('politicsuk_costs.mat')
% load('politicsuk_labels.mat')
% A = politicsuk_adjacency;
% C = politicsuk_costs;
% labels = politicsuk_labels';
% A = [A(1:2)];
% C = [C(1:2)];
addpath([pwd '/Datasets/3sources/'])  
load('sources_adj_density10.mat')
load('sources_costs_density10.mat')
load('sources_labels.mat')
A = sources_adjacencies;               
C = sources_weights;
labels = sources_labels;

m = numel(A);
k = max(labels);

% Run paramters
m_array = 1:m;          % no. of layers to run
num_runs = 15;           % no. of repetitions
ccr_array = zeros(numel(algs),numel(m_array),num_runs);
nmi_array = zeros(numel(algs),numel(m_array),num_runs);


% Data preprocessing
isol = [];                      % Find isolated nodes
for i = 1:m
        node_degrees = sum(A{i},2);
        isol = union(isol, find(node_degrees == 0));        
end

for i = 1:m                     % Remove isolated nodes
    A{i}(isol,:) = [];
    A{i}(:,isol) = [];
    C{i}(isol,:) = [];
    C{i}(:,isol) = [];
end
labels(isol) = [];
n = numel(labels);

for i = 1:m                     % Convert C to cost matrix from weight matrix
%     C{i} = 1./C{i};
%     C{i}(C{i} == inf) = 1e12;
    N = max(max(C{i}));
    C{i} = exp(-C{i}./N);
end

% Run comparison
for runs = 1:num_runs
    for alg_id = 1:numel(algs)
        for layers = 1:numel(m_array)
            % Run algorithms
            switch algs(alg_id)
                case 1
                disp('CRSP')
                [ccr_array(alg_id, layers, runs),nmi_array(alg_id, layers, runs), final_labels] = CRSP(A(1:layers),C(1:layers),labels', n, k, m_array(layers), b); 
                fprintf('CCR: %.2f\n', ccr_array(alg_id, layers, runs))
                fprintf('NMI: %.2f\n\n', nmi_array(alg_id, layers, runs))
                case 2
                disp('C-FE')
                [ccr_array(alg_id, layers, runs),nmi_array(alg_id, layers, runs), final_labels] = CFE(A(1:layers),C(1:layers),labels', n, k, m_array(layers), b); 
                fprintf('CCR: %.2f\n', ccr_array(alg_id, layers, runs))
                fprintf('NMI: %.2f\n\n', nmi_array(alg_id, layers, runs))
                case 3
                disp('SCML')
                addpath([pwd '/SCML/'])
                [ccr_array(alg_id, layers, runs), nmi_array(alg_id, layers, runs), final_labels] = SCML(A(1:layers), k, 0.5, labels);
                fprintf('CCR: %.2f\n', ccr_array(alg_id, layers, runs))
                fprintf('NMI: %.2f\n\n', nmi_array(alg_id, layers, runs))
                case 4
                disp('CPSC')
                sigma = zeros(1,m_array(layers));
                addpath([pwd '/coregularizedSC/'])  
                [ccr_array(alg_id, layers, runs), nmi_array(alg_id, layers, runs),final_labels] = spectral_pairwise_multiview(A(1:layers),m_array(layers),k,sigma,lambda_coreg, labels', num_iter);          
                fprintf('CCR: %.2f\n', ccr_array(alg_id, layers, runs))
                fprintf('NMI: %.2f\n\n', nmi_array(alg_id, layers, runs))
                case 5
                disp('CCSC')
                sigma = zeros(1,m_array(layers));
                addpath([pwd '/coregularizedSC/'])  
                [ccr_array(alg_id, layers, runs), nmi_array(alg_id, layers, runs),final_labels] = spectral_centroid_multiview(A(1:layers),m_array(layers),k,sigma,repmat(lambda_coreg,1,m), labels, num_iter);          
                fprintf('CCR: %.2f\n', ccr_array(alg_id, layers, runs))
                fprintf('NMI: %.2f\n\n', nmi_array(alg_id, layers, runs))
                case 6
                addpath([pwd '/MultiNMF/']) 
                disp('MultiNMF')
                options = [];
                options.maxIter = 75;
                options.error = 1e-6;
                options.nRepeat = 10;
                options.minIter = 30;
                options.meanFitRatio = 0.1;
                options.rounds = 30;
                options.WeightMode='Binary';
                options.varWeight = 0;
                options.kmeans = 1;
                options.Gaplpha=10;        % Change to 100?                      
                options.alpha=10;
                options.delta = 0.1;
                options.beta = 0;
                options.gamma = 2;
                options.K = k;
                options.alphas = ones(1,m);         % Equal weights to each layer
                [ccr_array(alg_id, layers, runs), nmi_array(alg_id, layers, runs)] = GMultiNMF(A(1:layers), k, C(1:layers), labels, options);
                fprintf('CCR: %.2f\n', ccr_array(alg_id, layers, runs))
                fprintf('NMI: %.2f\n\n', nmi_array(alg_id, layers, runs))
            end
        end
    end
end

if do_result_plot
    % Plot colors
    x =  [ 0    0.4470    0.7410;    0.8500    0.3250    0.0980;    0.9290    0.6940    0.1250; 0.4940    0.1840    0.5560; 0.4660    0.6740    0.1880;  0.3010    0.7450    0.9330;  0.6350    0.0780    0.1840]; 
    figure
    for i = 1:numel(algs)
        avg_ccr = mean(ccr_array(i,:,:),3);
        avg_nmi = mean(nmi_array(i,:,:),3);
        std_ccr = std(ccr_array(i,:,:),0,3);
        std_nmi = std(nmi_array(i,:,:),0,3);
        hold on;yyaxis left; errorbar(m_array, avg_ccr, std_ccr, '-', 'color', x(i,:),'DisplayName', sprintf('%s-CCR', alg_names(algs(i)))); 
        ylabel('CCR'); ylim([0,100])
        hold on;yyaxis right; errorbar(m_array, avg_nmi, std_nmi,'--', 'color', x(i,:),'DisplayName', sprintf('%s-CCR', alg_names(algs(i)))); 
        ylabel('NMI'); ylim([0,1])
    end
    title(sprintf('CRSP Comparisons: Nodes = %d, Clusters = %d, Layers = %d, Runs = %d',n,k,m, num_runs))
    xlabel('Layers'); 
    legend('Location','SouthEast')
 %saveas(gcf, [pwd '/figs/' sprintf('n%d_k%d_m%d_r%d.png', n,k,i,num_runs)])
end