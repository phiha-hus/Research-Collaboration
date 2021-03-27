
function [y] = starting1(funct,funct_d1,funct_d2,n,p,q,s,h,t0,y0,RTOL,ATOL)

     Y   = kron(ones(s,1),eye(n,n))*y0;
     dY  = kron(ones(s,1),eye(n,n))*y0;


     y(1:n,1)=y0;
     
     %c = 0:1/p:1;
     %[S] = startProc1(p,c,h);
     
     c = [1/4 0 -1/4];
     
     S = [ 1/4 0 0 1
          -1/4 1/4 0 1
           3/8 -7/8 1/4 1
           0 0 0 1
          1/2 2 -1/2 0
          8 -12 4 0];
  
     
     %internal stages -> Newton
     [Y,dY]= glm1_internal(funct,funct_d1,funct_d2,Y,dY,y,n,S,c,h,1,t0,p,q,s,1,RTOL,ATOL,1 );%% 1->2 

    % output quantities
     [y] =  glm1_external(funct,t0,y, dY,n,S,c,h,1,p,q,s,1 ); % 1->2

