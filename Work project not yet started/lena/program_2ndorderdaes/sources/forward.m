%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forward substitution b := L\b
%
% CALL  : [b]=forward(L,b,unit)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b=forward(L,b,unit)

if (nargin<3)
   unit=0;
end;

n=size(L,1);
if (unit==0)
   b(1,:)=b(1,:)/L(1,1);
end; % if
for i=2:n
    b(i,:)=b(i,:)-L(i,1:i-1)*b(1:i-1,:);
    if (unit==0)
       b(i,:)=b(i,:)/L(i,i);
    end; % if
end; % for
