% Main test file for C-RSP using external data (not SBM-generated)
% Anuththari Gamage
% 3/25/2018
clear;clc;close all

b = 0.02;                        % RSP parameter
do_result_plot = 1;              % To plot results

load('politicsuk_adjacency.mat')
load('politicsuk_costs.mat')
load('politicsuk_labels.mat')
A = politicsuk_adjacency;
C = politics_uk_costs;
labels = politicsuk_labels;

n = size(A{1},1);
m = numel(A);
k = 5;
ccr_array = zeros(1,m);
nmi_array = zeros(1,m);
m_array = 1:m;

for layers = 1:m
            % Run CRSP
            [ccr_array(i),nmi_array(i), final_labels] = CRSP(A,C,labels, n, k, m, b); 
            fprintf('CCR: %.2f\n', ccr_array(i))
            fprintf('NMI: %.2f\n\n', nmi_array(i))
end
    


if do_result_plot
 for i = 1:numel(m_array)
        figure;yyaxis left; plot(m_array, ccr_array); ylabel('CCR'); ylim([ 50,100])
        hold on;yyaxis right; plot(m_array, nmi_array); ylabel('NMI'); ylim([0,1])
        title(sprintf('CRSP: Nodes = %d, Clusters = %d, Layers = %d',n,k, m))
        xlabel('Layers'); 
        legend('CCR', 'NMI','Location','SouthEast')
 %saveas(gcf, [pwd '/figs/' sprintf('n%d_k%d_m%d_r%d.png', n,k,i,num_runs)])
 end
end
