function [G, P] = make_SBM(n,k,scaling_type,c,lambda)
% n is the number of nodes
% k is the number of communities
% scaling_type is either constant or logarithmic
% c is the constant fraction associated with the scaling
% lambda is the reduction percentage
% Brian Rappaport, 7/6/17


P = repmat(1:k,ceil(n/k),1);
P = P(1:n)';
if strcmp(scaling_type,'const')
    % odds if in same community is c/n, else is c(1-lambda)/n
    q = c*(1-lambda)/n;
elseif strcmp(scaling_type,'log')
    % odds if in same community is clog(n)/n, else is c(1-lambda)log(n)/n
    q = c*(1-lambda)*log(n)/n;
else
    error('scaling type must be ''const'' or ''log''');
end

% Build graph
indI = [];
indJ = [];
for i = 1:n
    I = rand(1,n);
    ind = P(i)==P;
    I(ind) = I(ind)*(1-lambda);
    f = find(I<q);
    indI(end+1:end+numel(f)) = f;
    indJ(end+1:end+numel(f)) = i;
end
G = sparse(indI,indJ,ones(numel(indI),1),n,n,numel(indI));
G = triu(G) + triu(G)';
end