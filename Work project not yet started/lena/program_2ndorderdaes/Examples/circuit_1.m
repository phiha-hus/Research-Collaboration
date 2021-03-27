
function f=circuit_1(t,v,dv)

R = 2;
L = 1/3;
C = 1.5;
vs=0;

f= [L*dv(4)-v(1)
    dv(2)-1/C*v(4)
    -R*v(4)+v(3)
    v(1)+v(2)+v(3)-vs];    