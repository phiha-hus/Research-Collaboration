%%========================================================================
% This is the testing example No 1
% Article Hoang - Hien - Phi 2020
% About triggering DAEs

clear all; close all; clc

E = [1 0 0; 0 1 0;0 0 0]; 

A = [-5 0 0; -2 -3 -2; -1 -3 -4] ; 

B = [0 0;1 0;0 1] ; 

alpha=0.5; sigma=0.5;

K = [2.0299 0.0400 0.0538 ; 
    1.3404 3.1254 0.2711];

x0 = [10 ; 2 ; 0.9801];

%=========================================================================
% Check strangeness-free 
syms t
[~,dim] = size(E) ; 
if rank(t*E-A)== dim
    sfree_flag = 1
else
    sfree_flag = 0
end
%=========================================================================

% Condition of the original system
t0 = 0; tf = 10; 
options = odeset('Mass',E,'RelTol',1e-4,'AbsTol',[1e-6 1e-10 1e-6]);
T_trig = [t0];

%=========================================================================
% Solving system on the time interval [0,10] and find the first triggering
% point

t_trig = t0 ; N = 1e+2; h = (tf-t_trig)/N; T = linspace(t_trig,tf,N);
x_trig = x0 ; 

% Update right hand side function
rhs_fun = @(t,x)A * x + B * K * x_trig;

[tspan,x] = ode15s(rhs_fun,T,x_trig,options);
x = x';

[trig_flag,t_trig,x_trig,idx] = find_trig_point(T,x,sigma);

T_plot = [];
X_plot = [];

if trig_flag == 1
    T_trig = [T_trig t_trig];
    
    T_plot = [T_plot T(1:idx)];
    X_plot = [X_plot x(:,1:idx)];
end

%=========================================================================
% Find the other triggering point
u_trig = K*x0;

while trig_flag == 1
   N = 1e+2; h = (tf-t_trig)/N; T = linspace(t_trig,tf,N); 
   
   % Update right hand side function
   rhs_fun = @(t,x)A * x + B * K * x_trig;
   
   [tspan,x] = ode15s(rhs_fun,T,x_trig,options);
   
   x = x';
   
   [trig_flag,t_trig,x_trig,idx] = find_trig_point(T,x,sigma) ;   
     
   if trig_flag == 1
    T_trig = [T_trig t_trig];
    
    % Update the piecewise solution on the interval [t_k,t_{k+1}]
    T_plot = [T_plot T(1:idx)];
    X_plot = [X_plot x(:,1:idx)];    
    
    u_trig = [u_trig K*x_trig];
   end      
end

plot_vidu
%=========================================================================






