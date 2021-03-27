function f=tanh_d2(t,x,dx,ddx,type)


if(strcmp(type,'bdf')) 

    f = [0 0
         0 0];
  
else
    f = [0 0
         0 0];
end