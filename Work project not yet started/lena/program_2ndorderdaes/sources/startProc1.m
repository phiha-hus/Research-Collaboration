% Starting procedure
%
% p order of method
% c vector of abscissae of lenght p+1
function [S] = startProc1(p,c)

if( length(c) ~= p+1)
    error('wrong input');
end

J = zeros(p+1,p+1);
C = zeros(p+1,p+1);
for i=1:p+1
    for j=1:p+1
        if( i == j+1)
            J(i,j)=1;
        end
        if j==1
            C(i,j)=1;
        else
            C(i,j)=c(i)^(j-1)/faculty(j-1);
        end
    end
end

B = J*inv(C);
A = C*J*inv(C);

U(:,1) = ones(p+1,1);

V = eye(p+1,1);
%V(1:2,1:2) = eye(2,2);

S = [ A U
      B V];
  
      