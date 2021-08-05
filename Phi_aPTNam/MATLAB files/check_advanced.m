% Check the advancedness and impulse-freeness of the system

function [is_non_advanced,is_impulse_free] = check_advanced(E,A,Ad)

is_non_advanced = 1; 
is_impulse_free = 1;

[n,m] = size(Ad);
k = floor(m/n) ;

Einf = [E' -A' zeros(n,n); zeros(n,n) E' -A'; zeros(n,2*n) E'; zeros(k*n,2*n) Ad'] ;
E1 = [E' -A'; zeros(n,n) E'; zeros(k*n,n) Ad'] ;
test = (rank(Einf)-rank(E1)-n == 0) ;

if test == 1
   disp('The system is non-advanced')
 else
   disp('The system is advanced')
   is_non_advanced = 0;
end

% Pre-processing to achieve prettier system
[U,S,V] = svd(E);
E = S ;
A = U' * A * V ;
% Ad = U' * Ad * V ;

E2 = [E(1:rank(E),:) ; A(rank(E)+1:end,:)] ;
if rank(E2)==n
    disp('System is impulse-free')
else
    disp('System is not impulse-free')
    is_impulse_free = 0;
end
