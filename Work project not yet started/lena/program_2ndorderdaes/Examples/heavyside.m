
function f=heavyside(t,x,dx,ddx)

n=3;

M=[1 1;t t];
C=[0 0;0 0];
K=[1 0;1+t 1];
g=[t^n*2*heaviside(t)-1;exp(t)];

f= M*ddx+C*dx+K*x-g;    