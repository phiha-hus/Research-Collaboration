function f=tanhrot_d3(t,x,dx,ddx,type)

R=[cos(t) sin(t);-sin(t) cos(t)];

if(strcmp(type,'bdf')) 

f = R*[ 1 1;t t];

else
    f = R*[ 1 1;t t];
    
end