clear;clc;close all

n = 100;
k = 2;
scaling_type = 'const';
c = 10;
lambda = 0.9;
m = 2;
 
[A, P] = make_mlSBM(n,k,m,scaling_type,c,lambda);

for i = 1:m
figure;spy(A(:,:,i))
title(sprintf('Avg.Density = %.2f',nnz(A(:,:,i))/(n*n)));
end

