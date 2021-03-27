%internal stages -> Newton
function     [Y,dY,converged] = glm1_internal(funct,funct_d1,funct_d2,dYn,y,n,M,c,h,j,t,p,s,k,RTOL,ATOL,start)

A = M(1:s,1:s);
U = M(1:s,s+1:s+k);

%starting value for Newton-Iteration
Y0 = dYn;

% maximal number of iterations
MAXIT = 20;

x(:,1)=Y0;   
dx(:,1)=Y0;
abbruch=1;

if( j>1)
    D = zeros(p+1,p+1);
    for jj=1:p+1
        D(jj,jj)=(h(j)/h(j-1))^(jj-1);
    end
    y=kron(D,eye(n,n))*y;
end

% Iteration counter
m=1; 

while(abbruch >= 0.33 && m<=MAXIT)
    
    dY = x(:,m);
    Y  = kron(U,eye(n,n))*y + h(j)*kron(A,eye(n,n))*dY;
    
    for i=1:s
       J(n*(i-1)+1:n*i,n*(i-1)+1:n*i) = feval(funct_d1,t+c(i)*h(j),Y(n*(i-1)+1:n*i),dY(n*(i-1)+1:n*i));
       K(n*(i-1)+1:n*i,n*(i-1)+1:n*i) = feval(funct_d2,t+c(i)*h(j),Y(n*(i-1)+1:n*i),dY(n*(i-1)+1:n*i));
    end

    for i=1:s
        F(n*(i-1)+1:n*i,1)  = feval( funct,t+c(i)*h(j),Y(n*(i-1)+1:n*i),dY(n*(i-1)+1:n*i));
    end
    
    DF = K + h(j)*kron(A,eye(n,n))*J; 
     
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

dY = x(:,m);


