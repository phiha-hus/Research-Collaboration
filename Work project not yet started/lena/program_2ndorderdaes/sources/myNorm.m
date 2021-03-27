%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted root mean square norm
%
% CALL  : [V]=myNorm(v,RTOL,ATOL)
%
% INPUT:    v       - vector v
%           RTOL    - relative tolerance
%           ATOL    - absolute tolerance
%
% OUTPUT:   V       - weighted norm of v
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [V]=myNorm(v,RTOL,ATOL)

n=length(v);

WT=zeros(n,1);
sum=0;
for i=1:n
    WT(i)=RTOL(i)*abs(v(i))+ATOL(i);
    sum=sum+(v(i)/WT(i))^2;
end

V =sqrt(1/n*sum);