%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation
%
% Interpolate values y(i-nn),..,y(i+1) 
%
% CALL  : [y]=interpol(t,y,m)
%
% INPUT:    t     - t=[t(1),t(2),...t(i),t(i+1)]
%           T     - time ti interpolate
%           y     - y=[y(1),y(2),..,y(i),y(i+1)]
%           m     - number of interpolated values
%
% OUTPUT:   y     - interpolated value at T
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function  z=interpol(t,y,T,m)

n=length(t);

sum=0;
for i=0:m
    prod=1;
    if(i>0)
        for j=0:i-1
            T1=T;
            T2=t(n-j);
            prod=prod*(T1-T2);
        end
    else
        prod=1;
    end
    
    sum = sum + prod*divDiff(t(n-i:n),y(:,n-i:n),i);
end

z=sum;