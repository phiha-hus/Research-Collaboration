% Supplementary file for article HuyPhi2021
% Test the exponential dichotomy of the electronic circuit with 
% Josephson junction - Ref. Ricardo Riaza [23]

R = 1; L = 1; I0 = 1; k = 1; G = 10;

% p = cos(k*phi2) as in the article
 p = -1 % generates exp. stab. sys
% p = 1 % generates exp. dicho. sys

E = [1 0 0 ;0 1 0 ;0 0 0]

A = [-R/L 0 1;0 0 1;1/L I0*k*p G]

eig(E,A)

syms t
det(t*E-A)

hA1 = A(1:2,1:2) - A(1:2,3)* inv(A(3,3)) * A(3,1:2)
hA3 = inv(A(3,3)) * A(3,1:2)    

% Check the phi Lipschitz of f_1
[eye(2)  -A(1:2,3)* inv(A(3,3))]