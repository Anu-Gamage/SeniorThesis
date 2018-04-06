% Create weight matrixes, adjacency matrices from 3 sources data using
% Pearson Coefficient (?) model
% Anuththari Gamage 4/6/2018

clear;clc;
load('sources.mat')
load('vocab.mat')

% Set up for affinity
n = 169;
load('overlap_id.mat') %load global article ids
probs = cell(1,3);      % construct term probability vectors for each source

for i = 1:3                 % for each source
   probs{i} = cell(1,n);        
   for j = 1:n              % for each article
        article_id = overlap(j);
        p = zeros(1,numel(vocab));
        terms =  sources{i}(sources{i}(:,2)==article_id,1);  %terms in article
        freq = sources{i}(sources{i}(:,2)==article_id,3); % frequency of term
        for k = 1:numel(terms)
           p(terms(k)) = freq(k);
        end
        p = p./sum(p);              % convert to probability distribution over all terms
        probs{i}{j} = p;
   end
end

% assign affinities
thresh = 0.0017;
sources_weights = cell(1,3);
sources_adjacencies = cell(1,3);
for source = 1:3 % for each source
    aff = zeros(n);
    for i = 1:n % loop over all entries in affinity matrix
        for j = 1:n
            aff(i,j) = probs{source}{i}*probs{source}{j}';       
        end
    end
    
    sources_adjacencies{source} = sparse(aff > thresh);  
   % spy(sources_adjacencies{source})
    aff(aff <= thresh) = 0;
    sources_weights{source} = sparse(aff);       %affinities
end
density = 100*(sum(nnz(sources_adjacencies{1}) + nnz(sources_adjacencies{2}) + nnz(sources_adjacencies{3}))/3)/(n*n);
save('sources_adj_pc.mat', 'sources_adjacencies')
save('sources_costs_pc.mat', 'sources_weights')
