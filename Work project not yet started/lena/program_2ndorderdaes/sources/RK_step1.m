%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge-Kutta step
%
% Performs one Runge-Kutta step
%
% CALL  :
% [t,h,y,y_2]=RK_step(funct,y,i,h,t,s,type,RTOL,ATOL)
%
% INPUT : funct - name of routine to determin F
%                 (Text-String)
%         y     - 
%         i     - number of current step
%         h     - current stepsize
%         t     - current time
%         s     - number of stages, s=2 !
%         RTOL  - relative error tolerance
%         ATOL  - absolute error tolerance
%         type  - type of Runge-Kutta method
%                 type = 'Gauss'
%                        'Radau'
%                        'Lobatto'
%
% OUTPUT: t    - time steps
%         h    - vector of stepsizes
%         y    - approximate solution at times t(i)
%         
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_new,converged]=RK_step1(funct,funct_d1,funct_d2,n,y,i,h,t,s,type,RTOL,ATOL)   

% starting value for Newton iteration
x0=zeros(s*n,1);

% maximal number of iterations
MAXIT = 20;

x(:,1)=x0;   
d(:,1)=x0;
abbruch=1;

% RK-coefficients
[A,b,c] = RK_coef1(s,type);

% Iteration counter
m=1;   
while(abbruch >= 0.33 && m<=MAXIT)
     
    dY = x(:,m);
   
    Y  = kron(ones(s,1),eye(n,n))*y(1:n,i) +h*kron(A,eye(n,n))*dY;
    
    for j=1:s
       F(n*(j-1)+1:n*j,1)  = feval( funct,t(i)+c(j)*h,Y(n*(j-1)+1:n*j),dY(n*(j-1)+1:n*j));
    
       K(n*(j-1)+1:n*j,n*(j-1)+1:n*j) = feval(funct_d1,t(i)+c(j)*h,Y(n*(j-1)+1:n*j),dY(n*(j-1)+1:n*j));
       M(n*(j-1)+1:n*j,n*(j-1)+1:n*j) = feval(funct_d2,t(i)+c(j)*h,Y(n*(j-1)+1:n*j),dY(n*(j-1)+1:n*j));
    end

    DF = M + h*kron(A,eye(n,n))*K;
    

    DF2=1/(h)*DF;
    F2=1/(h)*F;
    
    [LL,UU,PP]=lu(DF2);
    f_neu=-PP*F2;
    yy = forward(LL,f_neu);
    z1 = backward(UU,yy);
    
    d(:,m+1)=z1;
    x(:,m+1)=x(:,m)+d(:,m+1);
    
    m=m+1;
    if(m>2)
        rho = (myNorm(x(:,m)-x(:,m-1),RTOL,ATOL)/myNorm(x(:,2)-x(:,1),RTOL,ATOL))^(1/(m-1));
        abbruch = rho/(1-rho)*myNorm(x(:,m)-x(:,m-1),RTOL,ATOL);
    end
end % end while

if(m<MAXIT)
    converged = 0;
else
    converged = 1;
end

y_i = y(1:n,i);
y_ip1  = y_i + h*kron(b,eye(n,n))*x(:,m);

% the new approximation at t(i)+h
y_new(:,1)=y_ip1;


