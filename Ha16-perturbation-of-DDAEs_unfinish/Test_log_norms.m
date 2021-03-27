% Testing different concepts of the logarithmic norms
% Why the log norm mu[A^-1 * B] not= mu[B * A^-1], which is defined by 
% Higueras et al. in 1999

clear all; close all; clc

n = 2

for i=1:10

% Here we use A1 instead of A^-1 for convenience
A1 = diag(randn(n,1));

B = randn(n,n);

M1 = eye(n)+ A1 * B;

M2 = eye(n)+ B * A1;

defekt = norm(M1) - norm(M2)

check = (defekt > 0)

end