% Example taken from Xu.et.al'02
% To solve the following LMI
% E * P^T = P * E^T
% [A*P^T+P*A^T+Q   Ad * P^T
%   P * Ad^T       Q        ] < 0


clear all; close all; clc

E = [1 1 0; 1 -1 1; 2 0 1]
A = [1.5 0.5 1; -1 0 1; 0.5 0 1] 
Ad = [-1 0 -1; 1 -1 0.5; 0.3 0.5 -1]

A = A - 2 * E  % Shift the eigenvalue to have stability
% Theoretical A = A - E should work (eigenvalue method) but here LMI does
% not work at all

%% Check advancedness & impulse-freeness
[is_non_advanced,is_impulse_free] = check_advanced(E,A,Ad);
[E,A,Ad] = qz_transform(E,A,Ad)

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
