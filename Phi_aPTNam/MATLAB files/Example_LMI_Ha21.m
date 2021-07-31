% Example LMI
% Modified Example in Michiels'11
% serving the paper Ha21 on stability of Delay-DAEs
% @Copyright Phi Ha July, 2021

%%
clear all; close all; clc

E = [1 0 0; 0 0 1; 0 0 0] 
A = [0 -1/8 0; -1 1 0; 0 0 1] 

% Shift the eigenvalues - seems that Michiels is wrong
A = A - 2 * E

[EE,AA,Q,Z] = qz(E,A);
[EES,AAS,QS,ZS] = ordqz(EE,AA,Q,Z,'lhp');

EES
AAS
dEE = diag(EES)
dAA = diag(AAS)

a = 1/4;
% a = 3/4;

A_1d = [0 0 0; 0 -a 0; 0 0 0] 
A_2d = [0 0 0; 0 1/2 0; 0 0 0]

tau = [1 2]; h = [0 tau]  % Option 1
%tau = [0.99 2]; h = [0 tau]  % Option 2


%% Matrix Coefficients are motivated from 
% Cui.et.al TAC2017
% Index 2 system

clear all; close all; clc

alpha = 0;
%alpha = 0.1270; 

E = [-11.0000    1.0000     0 %0.7885
         0         0         alpha
         0         0         0]
  
A = [0.2000    0.6100    0.1891
   -1.0000    0.6000    0.5607
         0         0    0.2998]
     
Ad = [-1.0000   -0.2000   -1.5970
   -0.8000   -0.0100         0
         0         0         0]

syms s w
p1 = det(s*E-A-w*Ad);
% Very nice trick in MATLAB here. Computational error leads to p1 of degree
% 3, but in fact it is of degree 1.
disp('The matrix polynomial det(sE-A) reads')
vpaSols = vpa(p1,6)

[n,~] = size(E);

%% Test LMIs
setlmis([]) ;

structure_Q = [n,1]
Q = lmivar(1,structure_Q)

structure_P = [n,1]
P = lmivar(1,structure_P)

S = [null(E') zeros(n,rank(E))]

% LMI -1
lmiterm([-1,1,1,P],1,1)  % P > 0

% LMI -2
lmiterm([-2,1,1,Q],1,1)  % Q > 0


% LMI 1
lmiterm([1,1,1,P],A,E','s')
lmiterm([1,1,1,Q],A*S',eye(n),'s')
lmiterm([1,1,1,Q],1,1)

lmiterm([1,1,2,P],Ad,E')
lmiterm([1,1,2,Q],Ad * S',eye(n))

lmiterm([1,2,1,P],E,Ad')
lmiterm([1,2,1,Q],eye(n),S*Ad')

lmiterm([1,2,2,-Q],1,1)

% getlmis
LMISYS = getlmis ;

[t_min,X_feas] = feasp(LMISYS)

P = dec2mat(LMISYS,X_feas,P)
Q = dec2mat(LMISYS,X_feas,Q)

%% Verify the positive-definiteness of the matrices
    
% Cha hieu sao ket qua be the, toan phai scale
%P = 1e+11 * P
%Q = 1e+11 * Q

if min(eig(P))<0
    disp('P is not positive semi-definite')
elseif min(eig(Q))<0
    disp('Q is not positive semi-definite')
end

P2 = E * P + Q * S 

% 0 hieu sao lam lai kieu LMI thi no lai 0 ra
M = [A*P2+P2*A'+Q, Ad*P2 ; P2*Ad' , -Q]
%eig(M)

% Schur complement thi ra rat chuan
N = A*P2+P2*A' + Q + Ad*P2*inv(Q)*P2*Ad'
%eig(N)

if max(eig(M))>0
    disp('M is not negative definite')
elseif max(eig(N))>0
    disp('N (the Schur complement) is not negative definite')
end
   


