%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BDF method for the approximation of first order equations
%
% INPUT:  k     - order of the bdf method
%         T     - times T=[t(n),...,t(n-k)]
%         y     - y=[y(t_n),...,y(t_{n-k})]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_1] = bdf1(y,T,k)

sum=0;

for i=1:k
    prod = 1;
    if(i>1)
        for m=1:i-1
            prod = prod*(T(1)-T(m+1));
        end
    end   
    sum = sum + prod*divDiff(T,y,i);
end

y_1 = sum;