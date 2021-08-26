% testing the code TDS_STSBIL written by Wim Michiels
% serving the paper Ha21 on stability of Delay-DAEs
% Syntax must see in the function tds_create
% & compute_root_DDAE
% @Copyright Phi Ha July, 2021


%% Example by Phi - regular, impulse free
% Matrix Coefficients are taken from 
% Cui.et.al TAC2017
% Result is Unstable - So Cui is wrong?
% Only for sufficiently big tau (>=1.5) then stability holds.

clear all; close all; clc

E = [1 0;0 0]; 
A = [0.3640 0.3280;0  0]; 
Ad1 = [-0.4720   -0.0540;0 0]; 
Ad2 = [0 0;0 0];

% Update w.r.t parameters A^0_21 and A^1_21
alpha = -10 ; % parameter stands for A^0_21
beta = 0.01 ; % parameter stands for A^1_21

A(2,:) = alpha * A(1,:); A
Ad1(2,:) = alpha * Ad1(1,:) + beta * A(1,:); Ad1 
Ad2(2,:) = beta * Ad1(1,:); Ad2

s = 0
det(s * E - A - Ad1 - Ad2)

tau = [0.25 0.5]; h=[0 tau];
tds = tds_create({E},0,{A, Ad1, Ad2},h,'neutral')

%% 
% Reconstruct original Hessenberg DDAE
A_init = [A(1,:); alpha 0]
Ad_init = [Ad1(1,:); beta 0]

tau2 = [0.25]; h2 = [0 tau2];
tds2 = tds_create({E},0,{A_init, Ad_init},h2,'neutral')

s = 0
det(s * E - A_init - exp(-0.25 * s) * Ad_init)

%% Compute roots and plot
options=tdsrootsoptions;
% v = compute_roots_DDAE(tds2,options)
v = compute_roots_DDAE(tds,options)

vv = [v.l0; v.l1];
max(real(vv))

figure(1); clf;
plot(vv,'r*')
grid on
