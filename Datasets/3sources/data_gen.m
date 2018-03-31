clear;clc;close all

sources = cell(1,3);
ids = cell(1,3);

sources{1} = load('3sources_bbc.mtx');
sources{2} = load('3sources_guardian.mtx');
sources{3} = load('3sources_reuters.mtx');
ids{1} = load('3sources_bbc.docs');
ids{2} = load('3sources_guardian.docs');
ids{3} = load('3sources_reuters.docs');

for i = 1:3
    sources{i} = sources{i}(2:end,:); % trim metadata
    sources{i}(:,2) = ids{i}(sources{i}(:,2)); % Replace doc id with global doc ids
end

overlap = intersect(sources{1}(:,2), intersect(sources{2}(:,2), sources{3}(:,2)));  % 169 overlapping articles in all three sources
% Remove non-overlapping articles and sort by doc id
for i = 1:3
    sources{i} = sortrows(sources{i}(ismember(sources{i}(:,2),overlap),:),2);
end

thresh = 15;
n = numel(overlap);
sources_weights = cell(1,3);
sources_adjacencies = cell(1,3);
for source = 1:3 % for each source
    data = zeros(n);
    for i = 1:n % loop over all entries in adjacency matrix
        elems_i = sources{source}(sources{source}(:,2) == overlap(i),1);
        for j = 1:n
            % data(i,j) takes the number of overlapping words between articles i and j
            elems_j = sources{source}(sources{source}(:,2) == overlap(j),1);
            data(i,j) = numel(intersect(elems_i,elems_j));
        end
    end
    sources_adjacencies{source} = (data > thresh) + 0;
    sources_weights{source} = data;
    sources_weights{source}(data <= thresh ) = 0;
   
end
%%

topics = csvread('3sources.disjoint.clist',0,1);
topics = topics.*ismember(topics,overlap);
sources_labels = sort(topics(topics ~= 0));
for i = 1:numel(sources_labels)
    [r,~] = find(topics == sources_labels(i));
    sources_labels(i) = r;
end

save('sources_adj.mat', 'sources_adjacencies')
save('sources_costs.mat', 'sources_weights')
save('sources_labels.mat', 'sources_labels')