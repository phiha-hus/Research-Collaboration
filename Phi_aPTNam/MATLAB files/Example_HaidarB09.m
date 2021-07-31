% testing the code TDS_STSBIL written by Wim Michiels
% serving the paper Ha21 on stability of Delay-DAEs
% Syntax must see in the function tds_create
% & compute_root_DDAE
% @Copyright Phi Ha July, 2021


%% Example in Haidar/Boukas'09
clear all; close all; clc

E = [-1 2; -2 4]; 
A = [4.7 0.4; -4.9 0.8]; 
Ad1 = [0.7 -0.95; 1.1 -1.75]; 
Ad2 = [1 -0.8; 1.4 -1.3];

M = [1 0;-2 1]
E = M * E
A = M * A
Ad1 = M * Ad1
Ad2 = M * Ad2

tau = [0.2 2];

h=[0 tau];

% tds = tds_create({E},0,{A, Ad1, Ad2},h,'neutral')

%% Example by Phi - high index, non-advanced
% Motivated from Example in Haidar/Boukas'09
% Use this example to write article

%%%%%%%%%%%%%%%%%%%%
% Index 2 system

alpha = rand(2,1);
beta = rand(2,1);
gamma1 = rand(2,1);
gamma2 = rand(2,1);
epsilon = rand(1,1);

E = [E  alpha;zeros(1,3)]; 
A = [A  beta ;0 0 1];
Ad1 = [Ad1 gamma1; zeros(1,3)];
Ad2 = [Ad2 gamma2; zeros(1,3)];

MM = blkdiag(inv(M),epsilon)
E = MM * E
A = MM * A
Ad1 = MM * Ad1
Ad2 = MM * Ad2


E =[  -1.0000    2.0000    0.2648
   -2.0000    4.0000    0.8476
         0         0         0]


A = [4.7000    0.4000    0.1192
   -4.9000    0.8000    1.1783
         0         0    0.6473]


Ad1 = [0.7000   -0.9500    0.6456
    1.1000   -1.7500    1.7706
         0         0         0]


Ad2 = [1.0000   -0.8000    0.6393
    1.4000   -1.3000    1.8234
         0         0         0]
         
         
tds1 = tds_create({E},0,{A, Ad1, Ad2},h,'neutral')

%% Check index via QZ decomposition
[EE,AA,Q,Z] = qz(E,A);
[EES,AAS,QS,ZS] = ordqz(EE,AA,Q,Z,'lhp');

EES
AAS
AASd1 = QS * Ad1 * ZS
AASd2 = QS * Ad2 * ZS

dEE = diag(EES)
dAA = diag(AAS)

norm(EES - QS * E * ZS)
norm(AAS - QS * A * ZS)

tds2 = tds_create({EES},0,{AAS, AASd1, AASd2},h,'neutral')

%% Compute roots and plot
options=tdsrootsoptions;
v1 = compute_roots_DDAE(tds1,options)
vv1 = [v1.l0; v1.l1];
max(real(vv1))

v2 = compute_roots_DDAE(tds2,options)
vv2 = [v2.l0; v2.l1];
max(real(vv2))

figure(1); clf;
subplot(2,2,1)
plot(vv1,'r*')
grid on

subplot(2,2,2)
plot(vv2,'r*')
grid on

% print('Example13_Ha21','-depsc')
% print -depsc Example13_Ha21.eps
% ! epstopdf Example13_Ha21.eps


