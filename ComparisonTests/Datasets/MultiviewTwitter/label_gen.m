% Convert communities to label vector
% Anuththari Gamage 3/25/2018
clear;clc;close all

% load communities into 1x2 cell array before this (cell{1} = data, cell{2}
% = labels)
dataset = 'politicsuk';

load(sprintf('%s_communities.mat', dataset));
eval(sprintf('data = %s_communities;', dataset))
load(sprintf('%s_ids.mat', dataset));
eval(sprintf('id = %s_ids;', dataset));

labels = zeros(1,numel(id));
for i = 1:numel(id)
   x = id(i);
   j = 1;
   idxFound = 0;
   while j <= numel(data{1}) &&  idxFound == 0
         y = find(data{1}{j} == x);
         idxFound = ~isempty(y);  
         j = j+1;
   end
   labels(i) = j-1;
end
eval(sprintf('%s_labels = labels;', dataset));
save(sprintf('%s_labels.mat', dataset), sprintf('%s_labels', dataset))