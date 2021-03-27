function f=ex2Lin(t,x,dx,ddx)

M=[1 1;t t];
C=[0 0;0 0];
K=[1 0;1+t 1];
g=[sin(t)+t^2*cos(t);exp(t)];

f= M*ddx+C*dx+K*x-g;    