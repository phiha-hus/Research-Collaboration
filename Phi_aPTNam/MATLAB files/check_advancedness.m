% Checking the non-advancedness of the index 2 DDAE

E = [1 0;0 0]
A = [0 1;1 0];
[~,n] = size(E)

ep1 = 1
ep2 = 0

Ad = [0 ep1;0 ep2]

Einf = [E' -A' zeros(n,n); zeros(n,n) E' -A'; zeros(n,2*n) E'; zeros(n,2*n) Ad]

E1 = [E' -A'; zeros(n,n) E'; zeros(n,n) Ad']

test = (rank(Einf)-rank(E1)-n == 0)

if test == 1
   disp('The system is non-advanced')
 else
   disp('The system is advanced')
 end
 