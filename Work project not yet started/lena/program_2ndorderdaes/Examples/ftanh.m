function f=ftanh(t,x,dx,ddx)

n=3;
g=100;

M=[1 1;t t];
C=[0 0;0 0];
K=[1 0;1+t 1];
%g=[t^n*tanh(g*t+1)-1;exp(t)];
g=[0.5*t^2*(real(log(t))-1/2)-0.5*t^2;exp(t)];

f= M*ddx+C*dx+K*x-g;    