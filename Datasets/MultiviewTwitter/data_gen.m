% Convert data matrices to sparse .mat
% Anuththari Gamage 3/25/2018
clear;clc;close all

% Need to construct ids and communities first before this! 

dataset = 'politicsuk';
datatype = 'retweetedby';
data = load((sprintf('%s/%s-%s.mtx', dataset, dataset, datatype)));
data = data(2:end, :); % remove header

% replace ids with node index
load(sprintf('%s_ids.mat', dataset));
eval(['id = ',sprintf('%s_ids;', dataset)]);
rows = data(:,1);
cols = data(:,2);
for i = 1:numel(id)
    rows(rows==id(i)) = i;
    cols(cols==id(i)) = i;
end

data = sparse(rows, cols, data(:,3));
eval([sprintf('%s_%s = ',dataset, datatype),'data;']);
save(sprintf('%s_%s.mat',dataset, datatype), sprintf('%s_%s',dataset, datatype))