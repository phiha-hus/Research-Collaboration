%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% backward substitution b := U\b
%
% CALL  : [b]=backward(U,b,unit)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b=backward(U,b,unit)

if (nargin<3)
   unit=0;
end;

n=size(U,1);
if (unit==0)
   b(n,:)=b(n,:)/U(n,n);
end; % if
for i=n-1:-1:1
    b(i,:)=b(i,:)-U(i,i+1:n)*b(i+1:n,:);
    if (unit==0)
       b(i,:)=b(i,:)/U(i,i);
    end; % if
end; % for
