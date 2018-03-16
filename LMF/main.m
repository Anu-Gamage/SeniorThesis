% Main test file for LMF
% Anuththari Gamage
% 1/14/2017
clear;clc;close all
rng('shuffle')
n = 50;       % no. of nodes 
k = 6;        % no. of clusters
m = 1;        % no. of layers 
p = 3;        % SNR parameter (2-3 for strong signal, 1 for weak)
alpha = 1e-4; % Correction parameter
density = 0.0075; % percentage  (0.0075 for 1%)

params = 1:10;           % No. of layers
avgAcc = zeros(1, numel(params));
avgNmi = zeros(1, numel(params));
nodeDensity = zeros(1, numel(params));

for var = 1:numel(params)
    
    m = params(var);
    fprintf('Variable %d processing:\n', var)
    
    % Generate adjacency tensor (data)
    [A, labels] = mlsbm(n,k,m,p,density);

    % Average degree
    deg = zeros(1,m);
    for i = 1:m
        deg(i) = mean(sum(A(:,:,i),2));
    end
    nodeDensity(var) = 100*mean(deg/n);
    fprintf('Avg. Node Degree: %.2f\n', nodeDensity(var))

    % Run LMF

    numRuns = 50;
    acc = zeros(1, numRuns);
    nmi = zeros(1, numRuns);
    for i = 1:numRuns
        [acc(i),nmi(i), final_labels] = lmf(A, labels, k, alpha, 0);
    end
    % disp(acc)
    % disp(nmi)
    avgAcc(var) = sum(acc)/numRuns;
    avgNmi(var) = sum(nmi)/numRuns;
    fprintf('Avg. CCR: %.2f\n', avgAcc(var))
    fprintf('Avg. NMI: %.2f\n\n', avgNmi(var))

end
plot(params, avgAcc)
hold on; plot(params, 100*avgNmi)
avgDensity = sum(nodeDensity)/numel(params);
title(sprintf('Nodes = %d, Clusters = %d, Layers = %d, Avg. Degree = %.2f', n, k, m, avgDensity))
xlabel('Variable'); ylabel('Accuracy'); ylim([0,100])
legend('CCR', 'NMI')
saveas(gcf, sprintf('figs/n%d_k%d_m%d.jpg', n, k, m))
