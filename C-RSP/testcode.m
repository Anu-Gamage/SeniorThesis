clear;clc;close all
A = rand(1000);
n = 0;
 while ~((mean(sum(A,1)) == 1) & (mean(sum(A,2)') == 1))
        A = A./sum(A,1);
        A = A./sum(A,2);
        n = n+1
 end



% 
% % A =[1,2,3;4,5,6;7,8,9]
% % 
% % syms f(x)
% % f(x) = x^(1/3);
% % 
% % B = funm(A,f)
% 
% P(:,:,1) = [1,2,3;4,5,0;7,8,0];
% P(:,:,2) = 2*[1,2,0;0,5,6;7,8,9];
% P(:,:,3) = 3*[1,2,0;0,5,6;7,8,9];
% 
% B = combineP(P)
% 
% function newP = combineP(P)
%    n = size(P,1);
%    nz = sum(P ~= 0,3);
%    Pnz = P + (P==0);
%    
%    newP = ones(n);
%    for i = 1:size(P,3)
%     newP = newP.*Pnz(:,:,i);
%    end
%    
%    newP = arrayfun(@(P,nz) nthroot(P,nz), newP, nz);
% end
% 
% %    newP = zeros(n);
% %    for i = 1:n
% %        for j = 1:n
% %            if nz(i,j) ~= 0
% %                 x = 1;
% %                 for l = 1:m
% %                     if P(i,j,l) ~= 0
% %                        x = x*P(i,j,l);
% %                     end
% %                 end
% %                 newP(i,j) = x;
% %            end                      
% %        end
% %    end
% %    nz = nz + (nz == 0);         % To eliminate 0th roots
% %    newP = arrayfun(@(P,nz) nthroot(P,nz), newP, nz);