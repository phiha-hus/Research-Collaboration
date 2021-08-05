% Example LMIs taken from Michiels'11, Cui.et.al'17, BoukasHaidar'09,
% Xu.et.al'02
% serving the paper Ha21 on stability of Delay-DAEs
% @Copyright Phi Ha July, 2021

%% Matrix Coefficients are motivated from 
% Cui.et.al TAC2017
% Index 2 system

clear all; close all; clc

%alpha = 0;
alpha = 0.1270; 

E = [-11.0000    1.0000     0 %0.7885
         0         0         alpha
         0         0         0]
  
A = [0.2000    0.6100    0.1891
   -1.0000    0.6000    0.5607
         0         0    0.2998]
     
Ad = [-1.0000   -0.2000   -1.5970
   -0.8000   -0.0100         0
         0         0         0]

[n,~] = size(E);


%% Check advancedness & impulse-freeness
[is_non_advanced,is_impulse_free] = check_advanced(E,A,Ad);
[E,A,Ad] = qz_transform(E,A,Ad)

% % System obtained after applied QZ and throw something
E = [  -4.8020   -9.9469   -0.7885
         0    0.0000         0
         0         0         0]
A = [  0.6260   -0.1423   -0.1891
         0    1.1662    0.5607
         0         0    0.2998]
Ad = [ -0.6860   -0.7546    1.5970
    0.4202    0.6808   -0.0000
         0         0         0]

%% Test LMIs
[n,~] = size(E); 

setlmis([]) ;

structure_P = [n,1]
P = lmivar(1,structure_P)

structure_Q = [n,1]
Q = lmivar(1,structure_Q)

%R = lmivar(2,[n,n])   % R can be arbitrarily
structure_R = [n,1]
R = lmivar(1,structure_R)

S = [null(E) zeros(n,rank(E))]'
% E * S' check if E*S'==0
% P2 = E * P + R * S

% LMI -1
lmiterm([-1,1,1,P],1,1)  % P > 0

% LMI -2
lmiterm([-2,1,1,Q],1,1)  % Q > 0

%=========================================================================
% This is really the strange thing which I do not understand
% Without this inequality, the LMI returns P, Q but M is not negative
% definite
% LMI -3
lmiterm([-3,1,1,R],1,1)  % R > 0
%=========================================================================

% LMI 4
lmiterm([4,1,1,P],A,E','s')
lmiterm([4,1,1,R],A*S',1,'s')
lmiterm([4,1,1,Q],1,1)
%lmiterm([4,1,1,0],1)

lmiterm([4,1,2,P],Ad,E')
lmiterm([4,1,2,R],Ad * S',eye(n))

lmiterm([4,2,2,Q],-1,1)
%lmiterm([4,1,1,0],1)

% getlmis
LMISYS = getlmis ;

%options = [0,0,10,0,0];
%[t_min,X_feas] = feasp(LMISYS,options,0)
[t_min,X_feas] = feasp(LMISYS)

P = dec2mat(LMISYS,X_feas,P)
Q = dec2mat(LMISYS,X_feas,Q)
R = dec2mat(LMISYS,X_feas,R)

%% Verify the positive-definiteness of the matrices
    
if min(eig(P))<0
    disp('P is not positive semi-definite')
elseif min(eig(Q))<0
    disp('Q is not positive semi-definite')
end

P2 = E * P + R * S ;

norm(E * P2' - P2 * E')

% 0 hieu sao lam lai kieu LMI thi no lai 0 ra
M = [A*P2'+P2*A'+Q, Ad*P2' ; P2*Ad' , -Q]

% Schur complement thi ra rat chuan
N = A*P2'+P2*A' + Q + Ad*P2'*inv(Q)*P2*Ad'

if max(eig(M))>0
    disp('M is not negative definite')
    eig(M)
elseif max(eig(N))>0
    disp('N (the Schur complement) is not negative definite')
    eig(N)
end

