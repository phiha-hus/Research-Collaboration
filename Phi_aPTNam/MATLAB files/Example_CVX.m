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

E = [1 1 0; 1 -1 1; 2 0 1] ;
A = [1.5 0.5 1; -1 0 1; 0.5 0 1] ;
Ad = [-1 0 -1; 1 -1 0.5; 0.3 0.5 -1] ;
[n,~] = size(E) ;

% Pre-processing to achieve prettier system
[U,S,V] = svd(E);
E = S
A = U' * A * V
Ad = U' * Ad * V

cvx_begin sdp

variable P(n,n) symmetric
variable Q(n,n) symmetric
variable beta2 
variable gamma2

%vec(E * P' - P * E') >= zeros(n^2,1)
%vec(E * P' - P * E') <= zeros(n^2,1)

% [-beta2 * eye(n) E * P' - P * E';E * P' - P * E' -eye(n)] <= 0

% [A*P' + P*A' + Q, Ad * P'; P * Ad', -Q] <= 0
% minimize(beta2)

beta2 >= 1e-9
gamma2 >= 1e-9

P >= beta2 * eye(n)
[A*P + P*A + Q, Ad * P; P * Ad', -Q + gamma2 * eye(n)] <= 0

cvx_end
