% Main test file for RSP
% Anuththari Gamage
% 3/11/2018
clear;clc;close all

n = 500;       % no. of nodes 
k = 2;         % no. of clusters
b = 0.02;      % Tuning parameter for RSP
c = 5;         % Node degree
lambda = 0.9;
scaling_type = 'const';
doplot = 0;    % To plot data matrices


params = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0];               % Parameter being changed
acc = zeros(1, numel(params));
nmi = zeros(1, numel(params));

for var = 1:numel(params)
        
    c = params(var);              % Change parameter here
    fprintf('Variable %d processing:\n', var)
         
     % Generate adjacency tensor (data)
     [A, labels] = make_mlSBM(n,k,1,scaling_type, c, lambda);
     if doplot
       figure; spy(A);
     end
     
     %Remove isolated nodes
%     degs = sum(A,2);
%     x = (degs == 0);
%     A(x,:) = [];
%     A(:,x) = [];
%     n = size(A,1);
%     labels(degs == 0) = [];

    % Run RSP
    [acc(var), nmi(var), final_labels] = RSP(A, labels, n, k, b);
    fprintf('CCR: %.2f\n', acc(var))
    fprintf('NMI: %.2f\n\n', nmi(var))
end

figure;yyaxis left; plot(params, acc); ylabel('CCR'); ylim([0,100])
hold on;yyaxis right; plot(params, nmi); ylabel('NMI'); ylim([0,1])
title(sprintf('RSP: Nodes = %d, Clusters = %d', n, k))
xlabel('c'); 
legend('CCR', 'NMI','Location','SouthEast')
