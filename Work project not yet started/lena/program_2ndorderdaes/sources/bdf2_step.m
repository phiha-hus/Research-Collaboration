%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BDF step
%
% solves nonlinear system of equations with Newton's method
% DF*d=-F
%
% Iteration  x(m+1)=x(m)-c*inv(J)* F(x(m))
%
% CALL  :
% [y,x,converged,IER]=Newton(funct,t,h,x0,y,RTOL,ATOL,K,a)
%
% INPUT:    funct   - Name der aufzurufenden Routine, die die Matrizen M, D,
%                     S und Vektor f bestimmt
%
%           x0      - starting value for Newton-iteration
%           t       - current time
%           h       - current stepsize
%           K       - current order
%           a       - vector of bdf coefficients
%           RTOL    - relative tolerance
%           ATOL    - absolute tolerance
%
% OUTPUT: y         - approximation of solution at t(i+1)
%         x         - all iterates
%         coverged  - converged==0->Iteration has converged, converged==1->
%                     iteratiom fails to converge
%         IER       - IER=-1 ->DF singular, else regular
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function [ y,converged, IER ] = bdf2_step(funct,funct_d1,funct_d2,funct_d3,n,nL,KK,TT,t,h,i,y,RTOL,ATOL,exponent )

    % BDF coefficients of BDF(l,j,k)
    a = bdfk_coef(KK,TT);
    
    nq = n-nL;  % nL=1 !!

    %predictor value
    sum1=0;
    for ii=1:KK(1)+1
        sum1=sum1+phi_star(ii,i,y,t(1:i),h);
    end

    %starting value for Newton-Iteration
    x(:,1)=sum1;
    d(:,1)=sum1;

    % Newton-Iteration
    m=1;          
    maxIt = 4; 
    abbruch=1;
    converged = 0;

    while(abbruch >= 0.33 && m<=maxIt)
    
        y_0=y(1:n,:);
        y_1=y(n+1:n+nq,:);
    
        % Zusammenbau von F
        sum1 = 0;
        sum2 = 0;
        sum3 = 0;
        term_a = 0;
    
        for c=1:KK(1)
            sum1 = sum1 + a(c)*divDiff( t(i+1-c:i+1),[y_1(:,i+1-c:i),x(n+1:n+nq,m)],c);
            prod=1;
            for j=1:c-1
                prod = prod*(t(i+1)-t(i+1-j));
            end
            sum3 = sum3 + prod*divDiff(t(i+1-c:i+1),[y_0(1:nq,i+1-c:i),x(1:nq,m)],c);
            sum2 = sum2 + 1/(t(i+1)-t(i+1-c));
            prod2=1;
            for j=1:c
                prod2=prod2*(t(i+1)-t(i+1-j));
            end
            term_a=term_a+a(c)*1/prod2;
        end
        
        F =  [feval(funct,    t(i+1),x(1:n,m),x(n+1:n+nq,m), sum1 )
              x(n+1:n+nq,m)-sum3];
              
              
        D1=feval(funct_d1,t(i+1),x(1:n,m),x(n+1:n+nq,m), sum1);
        D2=feval(funct_d3,t(i+1),x(1:n,m),x(n+1:n+nq,m), sum1,'bdf')*term_a-feval(funct_d2,t(i+1),x(1:n,m),x(n+1:n+nq,m), sum1,'bdf');
        D3= eye(nq,nq+nL).*(-sum2);
        D4= eye(nq,nq);
        DF= [D1 D2
             D3 D4];

        if( det(DF) == 0) % J singular
            IER=-1;
        else
            IER=0; % J regular
        end; 
    
        % Gauss Zerlegung
        [L,U,P]=lu(DF);
        bb = P*(-F);
        %yy = L\bb;
        %d(:,m+1) = U\yy;
        yy = forward(L,bb);
        d(:,m+1) = backward(U,yy);
        
        
        x(:,m+1)=x(:,m)+d(:,m+1);
    
        m=m+1;
        if(m>2)
            rho = (myNorm(x(:,m)-x(:,m-1),RTOL,ATOL)/myNorm(x(:,2)-x(:,1),RTOL,ATOL))^(1/(m-1));
            abbruch = rho/(1-rho)*myNorm(x(:,m)-x(:,m-1),RTOL,ATOL);
            if(rho>0.9)
                converged = 1;
                y=x(:,m);
                return;
            end
        end
    end

    if(m > maxIt)
        converged = 1;
    end

    y=x(:,m);
    
    
    
    