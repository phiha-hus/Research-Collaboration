%% Example testing CVX
global n tau gamma

alpha = rand()
gamma = rand()

A = [-4 1;1 -2]
[n,~] = size(A);
cvx_begin sdp

variable P(n,n) symmetric
variable tau

P >= eye(n)
[A'*P + P*A + alpha*P + tau*gamma^2*eye(n), P;...
    P -tau*eye(n)] <= 0

cvx_end

cvx_status

P

tau 
%%  Example for my article
clear all; close all; clc

% E = [1 1 0; 1 -1 1; 2 0 1] ;
% A = [1.5 0.5 1; -1 0 1; 0.5 0 1] ;
% Ad = [-1 0 -1; 1 -1 0.5; 0.3 0.5 -1] ;

% Index 2 system
% alpha = 0;
 alpha = 0.1270; 

E = [-11.0000    1.0000    0.7885
         0         0         alpha
         0         0         0]  
A = [0.2000    0.6100    0.1891
   -1.0000    0.6000    0.5607
         0         0    0.2998]     
Ad = [-1.0000   -0.2000   -1.5970
   -0.8000   -0.0100         0
         0         0         0]
[n,~] = size(E) ;

% Pre-processing to achieve prettier system
% [U,S,V] = svd(E);
% E = S 
% A = U' * A * V 
% Ad = U' * Ad * V 

[E,A,Ad] = qz_transform(E,A,Ad);
E(2,3) = 0;

%
E
A
Ad
cvx_begin sdp

variable P(n,n) symmetric
variable Q(n,n) symmetric
variable beta2 
variable gamma2

S = [null(E) zeros(n,rank(E))]' 

beta2 >= 1e-9;
gamma2 >= 1e-9; 

Q >= beta2 * eye(n) 
P >= beta2 * eye(n)

[A*P*E' + E*P*A' + A*S'*Q + Q*S*A' + Q, Ad*P*E'+Ad*S'*Q; E*P*Ad'+Q*S*Ad', -Q + gamma2 * eye(n)] <= 0
cvx_end

%P = 1e+9 * P
P
eig(P)

%Q = 1e+9 * Q
Q
eig(Q)

M = [A*P*E' + E*P*A' + A*S'*Q + Q*S*A' + Q, Ad*P*E'+Ad*S'*Q; E*P*Ad'+Q*S*Ad', -Q + gamma2 * eye(n)]
eig(M)