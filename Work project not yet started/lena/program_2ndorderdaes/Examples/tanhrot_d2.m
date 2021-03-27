function f=tanhrot_d2(t,x,dx,ddx,type)

R=[cos(t) sin(t);-sin(t) cos(t)];

if(strcmp(type,'bdf')) 

    f = R*[0 0
         0 0];
  
else
    f = R*[0 0
         0 0];
end