%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation of Predictorpolynomial
%
%
% CALL  : p= phi_star(i,N,y,t,H)
%
%				
% OUTPUT: p    - value of predictor polynomial at t
%         
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p= phi_star(i,N,y,t,H)   % N=n


if(i==1)
    beta=1;
else
    p1=1;
    p2=1;
    for c=1:i-1
        p1=p1*psi(c,N+1,H);
        p2=p2*psi(c,N,H);
    end
   
    beta=p1/p2;
end


if(i==1)
    phi=y(:,N);
else
    p1=1;
    for c=1:i-1
        p1 = p1*psi(c,N,H);
    end
    
    phi=p1*divDiff(t,y,i-1-1);  
end


p = beta*phi;


