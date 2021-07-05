% Example 3 Phi_Nam_21
% Positive systems whose index is 2, but is still 
% positive and stable
% Matrix Coefficients are motivated from 
% Cui.et.al TAC2017

clear all; close all; clc

%%
% Generate partial monomila transformation
% Since the last coordinate x3 = 0, so it is partial
Q1 = [0 1.;1. 0; rand(1,2)];
Q = [Q1 rand(3,1)]
P = rand(3,3)

%%
E1 = [1.0 -11.0;0 0];
E = [E1 rand(2,1); zeros(1,3)];
rank(E);

A1 = [0.61 0.2;0.60 -1.0;0 0];
A = [A1 rand(3,1)];

Ad1 = [-0.2 -1.;-0.01 -0.80;0 0];
Ad = [Ad1 rand(3,1)];

%% 
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

% Check the non-impulse-free of the system
% pkg load symbolic # OCTAVE needs
syms s
s*E-A ;
p = det(s*E-A);
disp('The matrix polynomial det(sE-A) reads')
vpa(p,6)
% Function polynomialDegree Only occur in MATLAB 2021
% We can check by hand
% deg = polynomialDegree(p,s)  

%%
% Generate more beautiful system
% From the original system

E = P * E * Q 
A = P * A * Q 
Ad = P * Ad * Q 

syms s
s*E-A;
p1 = det(s*E-A);
% Very nice trick in MATLAB here. Computational error leads to p1 of degree
% 3, but in fact it is of degree 1.
disp('The matrix polynomial det(sE-A) reads')
vpaSols = vpa(p1,6)


%% Computing Drazin inverse using Jordan decomposition
% IS TERRIBLY BAD
alpha = 1 ;
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


