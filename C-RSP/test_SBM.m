clear;clc;close all

n = 500;
k = 2;
scaling_type = 'const';
%c = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0];
c = 5;
lambda = 0.9;

for var = 1:5
   
[G, P] = make_SBM(n,k,scaling_type,c,lambda);
figure;spy(G)


end