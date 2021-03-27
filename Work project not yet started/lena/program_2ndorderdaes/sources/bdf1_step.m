
function [x,converged,IER] = bdf1_step(funct,funct_d1,funct_d2,h,t,i,n,k,y,RTOL,ATOL)

    %predictor value
    x0=0;
    for ii=1:k+1
        x0=x0+phi_star(ii,i,y,t(1:i),h);
    end
    
    m=1;          
    maxIt = 20; 

    z(:,1)=x0;      % starting value
    d(:,1)=x0;
    abbruch=1;

    converged = 0;

    while(abbruch >= 0.33 && m<=maxIt)
        
       % Zusammenbau von F
        sum=0;
        sum2=0;
        for c=1:k
            prod=1;
            for j=1:c-1
                prod = prod*(t(i+1)-t(i+1-j));
            end
            sum = sum + prod*divDiff(t(i+1-c:i+1),[y(:,i+1-c:i),z(1:n,m)],c);
            sum2 = sum2 + 1/(t(i+1)-t(i+1-c));
        end

        F = feval(funct,t(i+1),z(1:n,m), sum );
        term=feval(funct_d2,t(i+1),z(1:n,m), sum).*sum2;
        DF= feval(funct_d1,t(i+1),z(1:n,m), sum) + term;
        
        
        if( det(DF) == 0) % J singular
            IER=-1;
        else
            IER=0; % J regular
        end; 
        
        % Gauss Zerlegung
        [L,U,P]=lu(DF);
        bb = P*(-F);
        yy = L\bb;
        d(:,m+1) = U\yy;
        
        z(:,m+1)=z(:,m)+d(:,m+1);
    
        m=m+1;
        if(m>2)
            rho = (myNorm(z(:,m)-z(:,m-1),RTOL,ATOL)/myNorm(z(:,2)-z(:,1),RTOL,ATOL))^(1/(m-1));
            abbruch = rho/(1-rho)*myNorm(z(:,m)-z(:,m-1),RTOL,ATOL);
            if(rho>0.9)
                converged = 1;
                x=z(:,m);
                return;
            end
        end
    end

    if(m > maxIt)
        converged = 1;
    end

    x = z(:,m);
    
    
    