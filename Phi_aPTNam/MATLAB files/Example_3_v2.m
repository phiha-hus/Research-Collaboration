% Example 3 Phi_Nam_21
% Positive systems whose index is 2, but is still 
% positive and stable
% Matrix Coefficients are motivated from 
% Cui.et.al TAC2017

clear all; close all; clc

%%
... Generate partial monomial transformation ...
... Since the last coordinate x3 = 0, so it is partial...
... Notice that Q must be a permutation matrix - in ...
... order to guarantee the positivity
    
Q1 = [0 1.;1. 0];
Q = blkdiag(Q1,rand(1,1))
%P = rand(3,3)
%Q = eye(3);
P = eye(3);

%% Create E, A and check their index, their impulse-freeness
E1 = [1.0 -11.0;0 0];
E = [E1 rand(2,1); zeros(1,3)];
rank(E);

A1 = [0.61 0.2;0.60 -1.0;0 0];
A = [A1 rand(3,1)];

n = size(E);
% Check the regularity of the system
for i = 1:(n+1)
    lambda = rand(1,1) ; 
    if rank(lambda * E - A) == n
        disp('Regular matrix pair')
        break
    else
        disp('Singular matrix pair')
    end
end

%%
% Generate more beautiful system
% From the original system

E = P * E * Q 
A = P * A * Q 


syms s
s*E-A;
p1 = det(s*E-A);
% Very nice trick in MATLAB here. Computational error leads to p1 of degree
% 3, but in fact it is of degree 1.
disp('The matrix polynomial det(sE-A) reads')
vpaSols = vpa(p1,6)


h = 1e-2; t0 = -200*h; tf = 200 * h;

for t = t0:h:tf
    
%Ad = [-0.2 -1. rand(1,1);-0.01 -0.80 rand(1,1);0 0 0];
Ad = [-0.2 -1. t;-0.01 -0.80 0;0 0 0];
Ad = P * Ad * Q 

%% Computing Drazin inverse using Jordan decomposition
% IS TERRIBLY BAD
alpha = 3;
U = alpha * E - A;  % VERY CAREFUL, DON'T USE INVERSE
hE = U\E
hA = U\A
hAd = U\Ad

rank(hE) - rank(E)

%n = size(E); 
%[V1,D] = eig(hE)
%k = 1  # In this example size of the dynamical part is 1

%hE_D = inv(V) * blkdiag(inv(D(1:k,1:k)),zeros(n-k)) * V 
%hE_D = inv(V) * blkdiag(1/D(1,1),zeros(2)) * V
%hE_D * hE

%% Computing using the algorithm

% pkg load symbolic # OCTAVE needs
hEd = Drazin_inverse(hE)
P = hEd * hE

Abar = hEd * hA
Ad_bar = hEd * hAd

hAD = Drazin_inverse(hA) % Drazin inverse of hA
K = ( P-eye(size(P)) ) * hAD * hAd 

if min(min(K>=0))==1
    t
    break
elseif t==tf
        disp('Change variable t')
end
  
end
    