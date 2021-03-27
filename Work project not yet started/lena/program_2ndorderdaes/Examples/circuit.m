
function f=circuit(t,v,dv,ddv)

R = 2;
L = 1/3;
C = 1.5;
vs=0;

f= [L*C*ddv(2)-v(1)
    -R*C*dv(2)+v(3)
    v(1)+v(2)+v(3)-vs];    
