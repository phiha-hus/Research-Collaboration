% This script test the controllability of 4 first order reformulation
% including companion form
% Example 2, Phi-Thuan'20

clear all; close all; clc

%M = [1 0;0 0]
%D = [1 0;0 0]
%K = [0 1;1 0]

M = zeros(2,2)
D = [0 1;0 0]
K = [0 0;1 0]
B = [0; 1]

E = [-K zeros(2,2);zeros(2,2) M]
A = [zeros(2,2) -K; -K -D]
%contrl(E,A)
eig(E,A)

B = [1 0]';
B1 = [zeros(2,1); B] ; 
rank([E B1])

T = null(E);
rank([E A*T B1])

disp('Test C-conty')
for i = 0:10
  rank([i*E-A B1])
end

