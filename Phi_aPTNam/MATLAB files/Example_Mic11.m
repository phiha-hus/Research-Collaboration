% testing the code TDS_STSBIL written by Wim Michiels
% serving the paper Ha21 on stability of Delay-DAEs
% Syntax must see in the function tds_create
% & compute_root_DDAE
% @Copyright Phi Ha July, 2021

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

tds1 = tds_create({E},0,{A, A_1d, A_2d},h,'neutral')

%% Modified Example in Michiels'11
% clear all; close all; clc

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

tds2 = tds_create({E},0,{A, A_1d, A_2d},h,'neutral')

%% Compute roots and plot
options=tdsrootsoptions;
v1 = compute_roots_DDAE(tds1,options)
vv1 = [v1.l0; v1.l1];
max(real(vv1))

v2 = compute_roots_DDAE(tds2,options)
vv2 = [v2.l0; v2.l1];
max(real(vv2))

figure(1); clf;
subplot(2,2,2)
plot(vv1,'r*')
grid on

subplot(2,2,1)
plot(vv2,'r*')
grid on

%print('Example14_Ha21','-depsc')
%print -depsc Example14_Ha21.eps
% ! epstopdf Example14_Ha21.eps



