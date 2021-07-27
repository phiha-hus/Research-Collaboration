% testing the code TDS_STSBIL written by Wim Michiels
% serving the paper on stability of Delay-DAEs that I'm writing
% Syntax must see in the function tds_create
% & compute_root_DDAE

clear all; close all; clc

%% Example in Haidar/Boukas'09
E = [-1 2; -2 4] 
A = [4.7 0.4; -4.9 0.8] 
A_1d = [0.7 -0.95; 1.1 -1.75] 
A_2d = [1 -0.8; 1.4 -1.3]

tau = [0.2 2];

h=[0 tau];

tds = tds_create({E},0,{A, A_1d, A_2d},h,'neutral')

%% Example by Phi - high index
% Constructed from Example_3_v2 @Copyright Phi Ha
% Matrix Coefficients are motivated from 
% Cui.et.al TAC2017
% Use this example to write article

clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%
% Index 2 system
% alpha = 0;
alpha = 0.1270; 

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


%%%%%%%%%%%%%%%%%%%%% 
% Original system, which is index 1
% E = [-11.0000  1.0000; 0 0]
% A = [0.2000    0.6100; -1.0000    0.6000]
% Ad = [-1.0000   -0.2000; -0.8000   -0.0100]

tau = 1; h = [0 tau];
tds = tds_create({E},0,{A, Ad},h,'neutral')

%% Check index via QZ decomposition
[EE,AA,Q,Z] = qz(E,A);
[EES,AAS,QS,ZS] = ordqz(EE,AA,Q,Z,'lhp');

EES
AAS
dEE = diag(EES)
dAA = diag(AAS)


%% Compute roots and plot
options=tdsrootsoptions;
v = compute_roots_DDAE(tds,options)

vv = [v.l0; v.l1];
max(real(vv))

figure(1); clf;
plot(vv,'r*')
grid on
