clear all; close all; clc

Q1 = [0 1.;1. 0];
Q = blkdiag(Q1,rand(1,1))
%P = rand(3,3)
%Q = eye(3);
P = eye(3);

%% Create E, A and check their index, their impulse-freeness
E1 = [1.0 -11.0;0 0];
E = [E1 rand(2,1); zeros(1,3)]
[n,~] = size(E)
rE = rank(E)

A1 = [0.61 0.2;0.60 -1.0;0 0];
A = [A1 rand(3,1)]

[EE,AA,Q,Z] = qz(E,A)

[EES,AAS,QS,ZS] = ordqz(EE,AA,Q,Z,'lhp')


%%
n = 5

Q = orth(rand(n,n))
P = orth(rand(n,n))

n2 = 4

dE = [rand(n2,1)- 2 * rand(n2,1) ; zeros(n-n2,1)];
dA = [rand(n2,1)- 2 * rand(n2,1) ; zeros(n-n2,1)];
E = diag(dE) ;
A = diag(dA) ;

E = Q * E * P
A = Q * A * P


[EE,AA,Q,Z] = qz(E,A);
[EES,AAS,QS,ZS] = ordqz(EE,AA,Q,Z,'lhp');

dEE = diag(EES)
dAA = diag(AAS)






