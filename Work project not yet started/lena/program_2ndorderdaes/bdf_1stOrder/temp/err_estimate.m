%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err] = err_estimate(c,y,h,p,s)

% A  = M(1:s,1:s);
% U  = M(1:s,s+2:s+k);
% bT = M(2*s+1,1:s );
% B  = M(2*s+2:2*s+s,1:s );
% vT = M(2*s+1,s+2:s+k );
% V  = M(2*s+2:2*s+s,s+2:s+k );

cp(:,1) = ones(s,1);
for i=1:p
    if( i>1)
        cp(:,i)=c.^(i-1)/faculty(i-1);
    end
    C(:,i)=cp(:,i);
    Ep1(i,1) = 1/faculty(p-i+1); 
end

%beta = inv(eye(n-1,n-1) )*(Ep1-B*c'.^p.\faculty(p));

%A = [ C 1; c^p/faculty(p) -beta];
%b = [0;1];


err = phiT*h*f + psiT*y;