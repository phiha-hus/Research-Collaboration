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

E = [1.0000  -11.0000; 0 0]
A = [0.61000  0.200; 0.60000  -1.0000]
Ad = [-0.2000 -1.0000; -0.0100  -0.8000]

tau = 1.5; h = [0 tau];
tds = tds_create({E},0,{A, Ad},h,'neutral')

%% Example by Phi - high index, non-advanced
% Constructed from Example_3_v2 
% Matrix Coefficients are motivated from 
% Cui.et.al TAC2017
% Use this example to write article

clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%
% Index 2 system
alpha = 0;
%alpha = 0.1270; 

E = [-11.0000    1.0000    0.7885
         0         0         alpha
         0         0         0]
  
A = [0.2000    0.6100    0.1891
   -1.0000    0.6000    0.5607
         0         0    0.2998]
     
Ad = [-1.0000   -0.2000   -1.5970
   -0.8000   -0.0100         0
         0         0         0]

syms s w
p1 = det(s*E-A-w*Ad);
% Very nice trick in MATLAB here. Computational error leads to p1 of degree
% 3, but in fact it is of degree 1.
disp('The matrix polynomial det(sE-A) reads')
vpaSols = vpa(p1,6)


tau = 1; h = [0 tau];
tds = tds_create({E},0,{A, Ad},h,'neutral')

%% Check index via QZ decomposition
[EE,AA,Q,Z] = qz(E,A);
[EES,AAS,QS,ZS] = ordqz(EE,AA,Q,Z,'lhp');

EES
AAS
AAdS = QS * Ad * ZS

dEE = diag(EES)
dAA = diag(AAS)

norm(EES - QS * E * ZS')
norm(AAS - QS * A * ZS')

%% Compute roots and plot
options=tdsrootsoptions;
v = compute_roots_DDAE(tds,options)

vv = [v.l0; v.l1];
max(real(vv))

figure(1); clf;
plot(vv,'r*')
grid on
