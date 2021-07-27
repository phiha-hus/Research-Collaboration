% testing the code TDS_STSBIL written by Wim Michiels
% serving the paper on stability of Delay-DAEs that I'm writing
% Syntax must see in the function tds_create
% & compute_root_DDAE

clear all; close all; clc

%% Example in Michiels'11
E = [1 0; 0 0] 
A = [0 -1/8; -1 1] 

% Shift the eigenvalues - seems that Michiels is wrong
A = A - 2 * E

a = 1/4;
% a = 3/4;

A_1d = [0 0; 0 -a] 
A_2d = [0 0;0 1/2]

tau = [1 2]; h = [0 tau]  % Option 1
%tau = [0.99 2]; h = [0 tau]  % Option 2

tds = tds_create({E},0,{A, A_1d, A_2d},h,'neutral')

%% Modified Example in Michiels'11

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

tds = tds_create({E},0,{A, A_1d, A_2d},h,'neutral')


%% Example by Phi - regular, impulse free
% Matrix Coefficients are taken from 
% Cui.et.al TAC2017
% Result is Unstable - So Cui is wrong?
% Only for sufficiently big tau (>=1.5) then stability holds.

% E = [1.0000  -11.0000; 0 0]
% A = [0.61000  0.200; 0.60000  -1.0000]
% Ad = [-0.2000 -1.0000; -0.0100  -0.8000]
% 
% tau = 1.5; h = [0 tau];
% tds = tds_create({E},0,{A, Ad},h,'neutral')

%% Compute roots and plot
options=tdsrootsoptions;

v = compute_roots_DDAE(tds,options)
vv = [v.l0; v.l1]
max(real(vv))

figure(1); clf;
plot(vv,'r*')
ylim([-10,10])
xlim([-4,4])
grid on


