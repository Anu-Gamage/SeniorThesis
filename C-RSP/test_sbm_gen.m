% Test SBM
% Anuththari Gamage
% 2/17/2018
clear;clc;close all

n = 10;
k = 2;
c = 5;
lambda = 0.9;
cin = c;
cout = ((1-lambda)*c);
seed = 10;

[A,conf_true] = sbm_gen(n,k,cin,cout,seed);
spy(A)