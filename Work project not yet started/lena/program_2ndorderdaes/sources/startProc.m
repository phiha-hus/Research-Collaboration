% Starting procedure
%
% p order of method
% c vector of abscissae of lenght p+1
function [S] = startProc(p,c,h)

% 
% if( length(c) ~= p+1)
%     error('wrong input');
% end
% 
% J = zeros(p+1,p+1);
% C = zeros(p+1,p+1);
% M = zeros(p+1,p+1);
% for i=1:p+1
%     for j=1:p+1
%         if( i == j+1)
%             J(i,j)=1;
%         end
%         if( i == j+2)
%             M(i,j)=1;
%         end
%         if j==1
%             C(i,j)=1;
%         else
%             C(i,j)=c(i)^(j-1)/faculty(j-1);
%         end
%     end
% end
% 
% At = C*M*inv(C);
% Bt = M*inv(C);
% A = C*J*inv(C);
% 
% Ut(:,1) = ones(p+1,1);
% Ut(:,2) = c';
% 
% Vt = zeros(p+1,2);
% Vt(1:2,1:2) = eye(2,2);
% 
% U = zeros(p+1,2);
% U(:,2) = 1/h*ones(p+1,1);
% 
% S = [ At Ut
%       A  U
%       Bt  Vt];


s = p+1;

%if( length(c) ~= s)
%    error('wrong input');
%end

J = zeros(s,s);
C = zeros(s,s);
M = zeros(s,s);
for i=1:s
    for j=1:s
        if( i == j+1)
            J(i,j)=1;
        end
        if( i == j+2)
            M(i,j)=1;
        end
        if j==1
            C(i,j)=1;
        else
            C(i,j)=c(i)^(j-1)/faculty(j-1);
        end
    end
end

At = C*M*inv(C);
Bt = M*inv(C);
A = C*J*inv(C);

Ut(:,1) = ones(s,1);
Ut(:,2) = c';

Vt = zeros(s,2);
Vt(1:2,1:2) = eye(2,2);

U = zeros(s,2);
U(:,2) = 1/h*ones(s,1);

S = [ At Ut
      A  U
      Bt  Vt];