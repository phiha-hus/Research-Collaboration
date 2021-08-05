% Example LMIs taken from Michiels'11, Cui.et.al'17, BoukasHaidar'09,
% Xu.et.al'02
% serving the paper Ha21 on stability of Delay-DAEs
% @Copyright Phi Ha July, 2021

%% Modified Example in Michiels'11
clear all; close all; clc

E = [1 0 0; 0 0 1; 0 0 0] 
A = [0 -1/8 0; -1 1 0; 0 0 1] 

% Shift the eigenvalues - seems that Michiels is wrong
A = A - 2 * E

[EE,AA,Q,Z] = qz(E,A);
[EES,AAS,QS,ZS] = ordqz(EE,AA,Q,Z,'lhp');

EES
AAS
dEE = diag(EES)
dAA = diag(AAS)

a = 1/4;
% a = 3/4;

A_1d = [0 0 0; 0 -a 0; 0 0 0] 
A_2d = [0 0 0; 0 1/2 0; 0 0 0]

tau = [1 2]; h = [0 tau]  % Option 1
%tau = [0.99 2]; h = [0 tau]  % Option 2
