%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  divided differences
%
% [z] = divDiff(t,y,i)
%
% INPUT:   t  -  t=[t_n,t_n-1,...,t_n-i]
%          y  -  y=[y_n,y_n-1,..,y_n-i]
%          i
%
% OUTPUT:  z  =  y[t_n,..,t_n-i]  divided differences
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z]=divDiff(t,y,i)

[n,m] = size(y);

if(n==0)
    z=0;
else


for k=1:i+1
    if( k==1 )
        for j=1:m
            dd((j-1)*n+1:j*n,1)=y(:,j);
        end
    else
        for j=1:m-k+1
            if(t(j+k-1)-t(j) == 0)
                dd((j-1)*n+1:j*n,k)=0;
            else
               dd((j-1)*n+1:j*n,k)=(dd(j*n+1:(j+1)*n,k-1)-dd((j-1)*n+1:j*n,k-1))/(t(j+k-1)-t(j));
            end
        end
    end
end

z = dd(1:n,i+1);

end