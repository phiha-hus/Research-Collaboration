
function [y] = starting(funct,funct_d1,funct_d2,funct_d3,n,p,q,s,h,t0,y0,y1,RTOL,ATOL)

     Y   = kron(ones(s,1),eye(n,n))*y0;
     dY  = kron(ones(s,1),eye(n,n))*y1;
     ddY = kron(ones(s,1),eye(n,n))*y1;
     

     y(1:n,1)=y0;
     y(n+1:2*n,1)=h*y1;
     
     c = 0:1/p:1;
     
     [S] = startProc(p,c,h);
     
     %internal stages -> Newton
     [Y,dY,ddY]= glm_internal(funct,funct_d1,funct_d2,funct_d3,Y,dY,ddY,y,n,S,c,h,1,t0,p,q,s,2,RTOL,ATOL,1 ); 

    % output quantities
     [y] =  glm_external(funct,t0,y, ddY,n,S,c,h,1,p,q,s,2 ); 
