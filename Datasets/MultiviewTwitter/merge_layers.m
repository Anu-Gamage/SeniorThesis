% Merge data layers into adjacency and cost tensors
% Anuththari Gamage 3/25/2018
clear;clc;close all

dataset = 'politicsuk';
layers = ["follows", "followedby", "retweets", "retweetedby", "mentions", "mentionedby"];

adj = cell(1,numel(layers));
cost = cell(1,numel(layers));
for i = 1:numel(layers)
    load(sprintf('%s_%s.mat', dataset, layers(i)));
    eval(sprintf('adj{%d} = (%s_%s~=0);', i, dataset, layers(i))) 
    eval(sprintf('cost{%d} = (%s_%s);', i, dataset, layers(i))) 
end

eval(sprintf('%s_adjacency = adj;', dataset))
save(sprintf('%s_adjacency.mat', dataset), sprintf('%s_adjacency', dataset))
eval(sprintf('%s_costs = cost;', dataset))
save(sprintf('%s_costs.mat', dataset), sprintf('%s_costs', dataset))