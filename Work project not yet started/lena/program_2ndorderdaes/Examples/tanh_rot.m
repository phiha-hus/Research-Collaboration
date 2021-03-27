function f=tanh_rot(t,x,dx,ddx)

n=3;
g=100;

R=[cos(t) sin(t);-sin(t) cos(t)];

M=[1 1;t t];
C=[0 0;0 0];
K=[1 0;1+t 1];
g=[t^n*tanh(g*t+1)-1;exp(t)];

f= R*M*ddx+R*C*dx+R*K*x-R*g;  