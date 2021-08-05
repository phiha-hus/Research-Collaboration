clear all; close all; clc

alpha = 0; 

E = [1   0    0
    0   1   alpha
    0         0    0];

Muoiad = [0 0 0
     0 0 0
     0 0 1];


A = [ -4.9816         0   0
         0    -5.2546    0
         0         0    -0.455]


Ad = [-1    0   0
   0   1    0
   0   0   -0.1291];
   
% Example taken from Xu.et.al'02
% To solve the following LMI
% E * P^T = P * E^T
% [A*P^T+P*A^T+Q+Muoiad*A+A'*Muoiad'   Ad * P^T+Muoiad*Ad
%   P * Ad^T       Q        ] < 0


% Q > 0
% P la ma tran bat ky

%% Test LMIs
[n,~] = size(E); 

setlmis([]) ;

structure_P = [n,1]
P = lmivar(1,structure_P)

structure_Q = [n,1]
Q = lmivar(1,structure_Q)

R = lmivar(2,[n,n])   % R can be arbitrarily

S = Muoiad;
% E * S' check if E*S'==0
% P2 = E * P + R * S

% LMI -1
lmiterm([-1,1,1,P],1,1)  % P > 0

% LMI -2
lmiterm([-2,1,1,Q],1,1)  % Q > 0

% LMI 1 PA+A'P'+Q+Muoiad*A+A'*Muoiad'
% = E'PA+A'P'E+R*S*A
lmiterm([1,1,1,P],E',A,'s')
lmiterm([1,1,1,R],1,S*A,'s')
lmiterm([1,1,1,Q],1,1)
lmiterm([1,1,1,0],Muoiad*A+A'*Muoiad')
lmiterm([1,1,1,1],1,eye(n))

lmiterm([1,1,2,P],E',Ad)
lmiterm([1,1,2,0],Muoiad*Ad)
lmiterm([1,1,2,R],1,S*Ad)

lmiterm([1,2,2,Q],-1,1)

% getlmis
LMISYS = getlmis ;

[t_min,X_feas] = feasp(LMISYS);

P = dec2mat(LMISYS,X_feas,P);
Q = dec2mat(LMISYS,X_feas,Q);
R = dec2mat(LMISYS,X_feas,R);

%% Verify the positive-definiteness of the matrices
    
% Cha hieu sao ket qua be the, toan phai scale
%P = 1e+9 * P
%Q = 1e+9 * Q

%if min(eig(P))<0
 %   disp('P is not positive semi-definite')
%elseif min(eig(Q))<0
  %  disp('Q is not positive semi-definite')
%end

P2 = E * P + R * S ; 

% 0 hieu sao lam lai kieu LMI thi no lai 0 ra
M = [A*P2'+P2*A'+Q, Ad*P2' ; P2*Ad' , -Q];
eig(M)

% Schur complement thi ra rat chuan
N = A*P2'+P2*A' + Q + Ad*P2'*inv(Q)*P2*Ad';
eig(N)

if max(eig(M))>0
   disp('M is not negative definite')
elseif max(eig(N))>0
   disp('N (the Schur complement) is not negative definite')
end
