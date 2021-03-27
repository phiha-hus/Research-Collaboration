%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error estmator by extrapolation
%
% interpolate values y(i-2),..,y(i) and extrapolate to i+1 
%
% CALL  : [z]=extrapol(t,y,m)
%
% INPUT:    t     - t=[t(1),t(2),...t(i),t(i+1)]
%
%           y     - y=[y(1),y(2),..,y(i),y(i+1)]
%           m     - number of interpolated values
%
% OUTPUT:   y     - extrapolated value at t(i+1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function  z=extrapol(t,y,m)

n=length(t);

sum=0;
for i=0:m
    prod=1;
    if(i>0)
        for j=0:i-1
            T1=t(n);
            T2=t(n-1-j);
            prod=prod*(T1-T2);
        end
    else
        prod=1;
    end
    
    sum = sum + prod*divDiff(t(n-1-i:n-1),y(:,n-1-i:n-1),i);
end

z=sum;
