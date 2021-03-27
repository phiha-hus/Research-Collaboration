function f=circuit_d2(t,i,di,ddi,type)

L = 1/3;
C = 1.5;
R = 2;

if(strcmp(type,'bdf')) 

    f = [ 0 0 0
        0 -R*C 0
        0 0 0];
  
else
    f = [ 0 0 0
      0 -R*C 0
      0 0 0];
end