%internal stages -> Newton
function     [Y,dY,ddY,converged] = glm_internal(funct,funct_d1,funct_d2,funct_d3,ddYn,y,n,M,c,h,j,t,s,k,RTOL,ATOL,start)

At = M(1:s,1:s);
Ut = M(1:s,s+1:s+k);
A  = M(s+1:s+s,1:s);
U  = M(s+1:s+s,s+1:s+k); 

%starting value for Newton-Iteration
Y0 = ddYn;

% maximal number of iterations
MAXIT = 20;

x(:,1)=Y0;   
dx(:,1)=Y0;
abbruch=1;

if( j>1)
    D = zeros(k,k);
    for jj=1:k
        D(jj,jj)=(h(j)/h(j-1))^(jj-1);
    end
    y=kron(D,eye(n,n))*y;
end

% Iteration counter
m=1; 

while(abbruch >= 0.33 && m<=MAXIT)
    
    ddY = x(:,m);
    
    Y  = kron(Ut,eye(n,n))*y + h(j)^2*kron(At,eye(n,n))*ddY;
    dY = kron(U,eye(n,n))*y  + h(j)*kron(A,eye(n,n))*ddY;   
    
    for i=1:s
       J(n*(i-1)+1:n*i,n*(i-1)+1:n*i) = feval(funct_d1,t+c(i)*h(j),Y(n*(i-1)+1:n*i),dY(n*(i-1)+1:n*i),ddY(n*(i-1)+1:n*i));
       K(n*(i-1)+1:n*i,n*(i-1)+1:n*i) = feval(funct_d2,t+c(i)*h(j),Y(n*(i-1)+1:n*i),dY(n*(i-1)+1:n*i),ddY(n*(i-1)+1:n*i),'glm');
       L(n*(i-1)+1:n*i,n*(i-1)+1:n*i) = feval(funct_d3,t+c(i)*h(j),Y(n*(i-1)+1:n*i),dY(n*(i-1)+1:n*i),ddY(n*(i-1)+1:n*i),'glm');
    end

    for i=1:s
        F(n*(i-1)+1:n*i,1)  = feval( funct,t+c(i)*h(j),Y(n*(i-1)+1:n*i),dY(n*(i-1)+1:n*i),ddY(n*(i-1)+1:n*i));
    end
    
    DF = L + h(j)^2*kron(At,eye(n,n))*J + h(j)*kron(A,eye(n,n))*K; 
    
    DF=1/(h(1)^2)*DF;
    F=1/(h(1)^2)*F;
     
    [LL,UU,PP]=lu(DF);
    f_neu=-PP*F;
    yy = forward(LL,f_neu);
    if(start == 1)
        z1 = backward(UU,yy,1);
    else
        z1 = backward(UU,yy);
    end
    
    dx(:,m+1)=z1;
    x(:,m+1)=x(:,m)+dx(:,m+1);
    
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

ddY = x(:,m);

