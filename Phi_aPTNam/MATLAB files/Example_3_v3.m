% Example 3 Phi_Nam_21
% Positive systems whose index is 2, but is still 
% positive and stable
% Matrix Coefficients are motivated from 
% Cui.et.al TAC2017

clear all; close all; clc

E = [-11.0000    1.0000    0.1521
         0         0    0.9365
         0         0         0]
A = [0.2000    0.6100    0.9236
   -1.0000    0.6000    0.4683
         0         0    0.7722]
Ad = [-1.0000   -0.2000   -1.9298
   -0.8000   -0.0100         0
         0         0         0]

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

syms s
s*E-A;
p1 = det(s*E-A);
% Very nice trick in MATLAB here. Computational error leads to p1 of degree
% 3, but in fact it is of degree 1.
disp('The matrix polynomial det(sE-A) reads')
vpaSols = vpa(p1,6)


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

[Up,Sp,Vp] = svd(P);
Sp

A1 = Up' * Abar * Vp

% H1 = Up' * H * Up has the same stability property
... as H
    
H1 = blkdiag(A1(1,1)/Sp(1,1),diag(-rand(2,1))) ;

H = Up * H1 * Up'
eig(H)  %-0.5499    -0.3264    -0.1304

n = 3;
I = eye(n);

H_bar = [Ad_bar+H Ad_bar; K K-I] ;
eig(H_bar)

matrices_in_Example3

    