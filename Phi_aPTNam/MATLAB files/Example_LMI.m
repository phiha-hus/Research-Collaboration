% Example taken from Xu.et.al'02
clear all; close all; clc

E = [1 1 0; 1 -1 1; 2 0 1]
A = [1.5 0.5 1; -1 0 1; 0.5 0 1]
Ad = [-1 0 -1; 1 -1 0.5; 0.3 0.5 -1]

%% Check the advancedness and impulse-freeness of the system
[n,~] = size(E);

Einf = [E' -A' zeros(n,n); zeros(n,n) E' -A'; zeros(n,2*n) E'; zeros(n,2*n) Ad] ;
E1 = [E' -A'; zeros(n,n) E'; zeros(n,n) Ad'] ;
test = (rank(Einf)-rank(E1)-n == 0) ;

if test == 1
   disp('The system is non-advanced')
 else
   disp('The system is advanced')
end

% Pre-processing to achieve prettier system
[U,S,V] = svd(E);
E = S 
A = U' * A * V
Ad = U' * Ad * V

E2 = [E(1:rank(E),:) ; A(rank(E)+1:end,:)] ;
if rank(E2)==n
    disp('System is impulse-free')
else
    disp('System is not impulse-free')
end

%% Test LMIs
setlmis([]) ;

structure_Q = [n,1]
Q = lmivar(1,structure_Q)

structure_P = [1,1]
P = lmivar(1,structure_P)

% LMI 1
lmiterm([-1,1,1,Q],1,1)  % Q>0
lmiterm([-1,1,1,P],E,eye(n))  % E * P>0

% LMI 2
lmiterm([2,1,1,P],A,eye(n),'s')
lmiterm([2,1,1,Q],1,1)

lmiterm([2,1,2,P],1,Ad')
lmiterm([2,2,1,P'],Ad,eye(n))

lmiterm([2,2,2,-Q],1,1)

% getlmis
LMISYS = getlmis ;

[t_min,X_feas] = feasp(LMISYS)

P = dec2mat(LMISYS,X_feas,P)
Q = dec2mat(LMISYS,X_feas,Q)

% Cha hieu sao ket qua be the, toan phai scale
P = 1e+11 * P
Q = 1e+11 * Q
eig(P*E)

% 0 hieu sao lam lai kieu LMI thi no lai 0 ra
M = [A*P+P*A'+Q, Ad*P ; P*Ad' , -Q]
eig(M)

% Schur complement thi ra rat chuan
N = A*P+P*A'+Q + Ad*P * inv(Q) * P*Ad'
eig(N)

