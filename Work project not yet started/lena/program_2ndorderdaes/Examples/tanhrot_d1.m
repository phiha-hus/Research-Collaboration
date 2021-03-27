function f=tanhrot_d1(t,x,dx,ddx) % Ableitung nach x

R=[cos(t) sin(t);-sin(t) cos(t)];

f = R*[1 0; 1+t 1];