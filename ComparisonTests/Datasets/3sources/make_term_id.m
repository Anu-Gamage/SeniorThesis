% Generates global ids for all terms in articles in the three sources
% instead of local ids
% Anuththari Gamage

clear;clc;close all
addpath([pwd '\original_data\'])
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

load('terms_id.mat')            % load global id of terms for each source
load('vocab.mat')               % load global vocabulary

% Replace term id in sources with global term ids
for i = 1:3
    for j = 1:size(sources{i},1)  
    id = sources{i}(j,1);
    sources{i}(j,1) = terms_id{i}(id);
    end
end
save('sources.mat', 'sources')      % all data with global ids
save('overlap_id.mat', 'overlap')

clear;clc
load('terms.mat')
load('vocab.mat')
terms_id = cell(1,3);

for j = 1:3
new_terms = terms{j};
id = zeros(1,numel(new_terms));

for i = 1:numel(new_terms)
    id(i) = find(vocab == new_terms(i));
end
terms_id{j} = id;
end

save('terms_id.mat', 'terms_id')