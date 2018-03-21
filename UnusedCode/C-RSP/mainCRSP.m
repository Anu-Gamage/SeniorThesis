% Main test file for C-RSP
% Anuththari Gamage
% 3/11/2018
clear;clc;%close all

n = 100;                        % no. of nodes 
k = 2;                          % no. of clusters
m_array = [2,3];                % no. of layers
b = 0.02;                       % Tuning parameter for RSP
c = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0,15.0, 20.0];               % Varying node degree
lambda = 0.9;
scaling_type = 'const';
doplot = 0;                     % To plot data matrices
num_runs = 5;
ccr_array = zeros(num_runs, numel(c), numel(m_array));
nmi_array = zeros(num_runs, numel(c), numel(m_array));


for runs = 1:num_runs
    % Generate adjacency tensor (data)
    data = zeros(n,n,max(m_array),numel(c));
    
    for i = 1:numel(c)
        [data(:,:,:,i), labels] = make_mlSBM(n,k,max(m_array),scaling_type, c(i), lambda);
    end

    for layers = 1:numel(m_array)        % Varying no. of layers
    m = m_array(layers);

    acc = zeros(1, numel(c));
    nmi = zeros(1, numel(c));

        for degree = 1:numel(c)        % Varying node degree

            A = data(:,:,1:m,degree);   % Select relevant tensor
            fprintf('Variable %d processing:\n', degree)

             if doplot
               figure; 
               for i = 1:m
                  subplot(1,m,i);spy(A(:,:,i)); title(sprintf('Layer %d', i))
               end
             end

            % Run CRSP
            [acc(degree), nmi(degree), final_labels] = CRSP(A,A,labels', n, k,m, b); % Cost matrix = A
            fprintf('CCR: %.2f\n', acc(degree))
            fprintf('NMI: %.2f\n\n', nmi(degree))
        end
        
        ccr_array(runs,:, layers) = acc;
        nmi_array(runs,:, layers) = nmi;
    end
end

 for i = 1:numel(m_array)
        avg_ccr = mean(ccr_array(:,:,i));
        avg_nmi = mean(nmi_array(:,:,i));
        std_ccr = std(ccr_array(:,:,i));
        std_nmi = std(nmi_array(:,:,i));
        figure;yyaxis left; errorbar(c, avg_ccr, std_ccr); ylabel('CCR'); ylim([ 0,100])
        hold on;yyaxis right; errorbar(c, avg_nmi, std_nmi); ylabel('NMI'); ylim([0,1])
        title(sprintf('CRSP: Nodes = %d, Clusters = %d, Layers = %d', n, k,m_array(i)))
        xlabel('c'); 
        legend('CCR', 'NMI','Location','SouthEast')
end