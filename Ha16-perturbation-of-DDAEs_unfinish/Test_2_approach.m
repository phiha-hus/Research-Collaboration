% Test 2 approaches 
% for decoupling the system for stability investigation
% of linear DDAEs


n=3;

M = randn(n,n);
L = randn(n,n);
N = randn(n,n);


E = [eye(3) N; zeros(3,2*n)]

[U,S,V] = svd(E)






























