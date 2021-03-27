%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM-Verfahren, variable-step variable-order GLM
%
% Loesen der differentiell-algebraischen Gleichung 2. Ordung der Form
% M y''(t)=D y'(t) + S y(t) + f(t,y,y')
%
% mit GLM der Ordnung Ordnung k
%
% CALL  : [t,h,y,y_est,y_interpol]=bdfk2dae(funct,y0,y1,t0,a,N,K,RTOL,ATOL,H0)
%
% INPUT : funct - Name der aufzurufenden Routine, die die Matrizen M, D, S und Vektor f bestimmt
%                 (Text-String)
%         y0,y1 - initial value
%         t0    - initial time         t\in[t0,t0+a]
%         a     - length of interval   t\in[t0,t0+a]
%         N     - number of steps
%         K     - maximal order
%
% OUTPUT: t    - time steps t=[t0,t0+h,...,t0+(N-1)*h,t0+N*h=t0+a]
%         h    - stepsizes
%         y    - approximation of solution at t(i)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,h,y]=glm1dae_test(funct,funct_d1,funct_d2,y0,t0,a,p,q,k,s,RTOL,ATOL,h)

t(1)=t0;
TOUT=t0+a;
[n,m1]=size(y0);

% % Starting method
%  YS   = zeros(s*n,1);% kron(ones(s,1),eye(n,n))*y0;
%  dYS  = zeros(s*n,1);%kron(ones(s,1),eye(n,n))*y0;
%  yS(1:n,1)=y0;
%          
%  c = [1/4 0 -1/4];
%      
%  S = [ 1/4  0   0   1
%       -1/4  1/4 0   1
%        3/8 -7/8 1/4 1
%        0    0   0   1
%        1/2  2  -1/2 0
%        8   -12  4   0];
%        
% AS = S(1:s,1:s);
% US = S(1:s,s+1:s+1);
% BS = S(s+1:s+s,1:s );
% VS = S(s+1:s+s,s+1:s+1 );
% 
% %starting value for Newton-Iteration
% Y0 = dYS;
% 
% MAXIT = 6;
% x(:,1)=Y0;   
% dx(:,1)=Y0;
% abbruch=1;
% m=1; 
% while(abbruch >= 0.33 & m<=MAXIT)
%     
%     dYS = x(:,m);
%     YS  = kron(US,eye(n,n))*yS + h*kron(AS,eye(n,n))*dYS;
%     
%     for j=1:s
%        J(n*(j-1)+1:n*j,n*(j-1)+1:n*j) = feval(funct_d1,t0+c(j)*h,YS(n*(j-1)+1:n*j),dYS(n*(j-1)+1:n*j));
%        K(n*(j-1)+1:n*j,n*(j-1)+1:n*j) = feval(funct_d2,t0+c(j)*h,YS(n*(j-1)+1:n*j),dYS(n*(j-1)+1:n*j));
%        F(n*(j-1)+1:n*j,1)  = feval( funct,t0+c(j)*h,YS(n*(j-1)+1:n*j),dYS(n*(j-1)+1:n*j));
%     end
%     
%     DF = K + h*kron(AS,eye(n,n))*J; 
%      
%     [LL,UU,PP]=lu(DF);
%     f_neu=-PP*F;
%     yy = forward(LL,f_neu);
%     %dx(:,m+1) = backward(UU,yy,1);
%     dx(:,m+1) = backward(UU,yy);
%  
%     x(:,m+1)=x(:,m)+dx(:,m+1);
%     
%     m=m+1;
%     if(m>2)
%         rho = (myNorm(x(:,m)-x(:,m-1),RTOL,ATOL)/myNorm(x(:,2)-x(:,1),RTOL,ATOL))^(1/(m-1));
%         abbruch = rho/(1-rho)*myNorm(x(:,m)-x(:,m-1),RTOL,ATOL);
%     end
% end % end while
% 
% dYS = x(:,m);
% 
% % output quantities
% y(:,1) = kron(VS,eye(n,n))*yS + h*kron(BS,eye(n,n))*dYS;
     

y(:,1)=[1;-1;2;-10*h(1);10*h(1);-20*h(1);100*h(1)^2;-100*h(1)^2;200*h(1)^2];

Y(:,1)    = kron(ones(s,1),eye(n,n))*y0;
dY(:,1)   = kron(ones(s,1),eye(n,n))*y0;

i=1;
while(i<700 && t(i)<=TOUT)
    t(i+1) = t(i)+h;

    [M,c] = glm1coef(p,q,k,s);
    A = M(1:s,1:s);
    U = M(1:s,s+1:s+k);
    B = M(s+1:s+s,1:s );
    V = M(s+1:s+s,s+1:s+k );
           
    %starting value for Newton-Iteration
    Y0 = dY(:,i);   % zeros

    % maximal number of iterations
    MAXIT = 6;

    x(:,1)=Y0;   
    dx(:,1)=Y0;
    abbruch=1;
    m=1; 
    while(abbruch >= 0.33 && m<=MAXIT)
    
        dYn = x(:,m);
    
        Yn = kron(U,eye(n,n))*y(:,i) + h*kron(A,eye(n,n))*dYn;
    
        for j=1:s
            J(n*(j-1)+1:n*j,n*(j-1)+1:n*j) = feval(funct_d1,t(i)+c(j)*h,Yn(n*(j-1)+1:n*j),dYn(n*(j-1)+1:n*j));
            K(n*(j-1)+1:n*j,n*(j-1)+1:n*j) = feval(funct_d2,t(i)+c(j)*h,Yn(n*(j-1)+1:n*j),dYn(n*(j-1)+1:n*j));
            F(n*(j-1)+1:n*j,1)  = feval( funct,t(i)+c(j)*h,Yn(n*(j-1)+1:n*j),dYn(n*(j-1)+1:n*j));
        end
    
        DF = K + h*kron(A,eye(n,n))*J; 
     
        [LL,UU,PP]=lu(DF);
        f_neu=-PP*F;
        yy = forward(LL,f_neu);
        dx(:,m+1) = backward(UU,yy);
        
        x(:,m+1)=x(:,m)+dx(:,m+1);
    
        m=m+1;
        if(m>2)
            rho = (myNorm(x(:,m)-x(:,m-1),RTOL,ATOL)/myNorm(x(:,2)-x(:,1),RTOL,ATOL))^(1/(m-1));
            abbruch = rho/(1-rho)*myNorm(x(:,m)-x(:,m-1),RTOL,ATOL);
        end
    end % end while

    dY(:,i+1) = x(:,m);
    Y(:,i+1) = Yn;

    y(:,i+1) = kron(V,eye(n,n))*y(:,i) + h*kron(B,eye(n,n))*dY(:,i+1);

    i=i+1;
end


















