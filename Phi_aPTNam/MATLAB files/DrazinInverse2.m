function DrazinInverse2a = DrazinInverse2(a)
% Compute the Drazin Inverse of a matrix 0a0 using the limited
% algorithm.
%Need not computing the index of 0a0.
global q1 q2 s1 s2
logr = 1; logs = 1;
[m,n] = size(a);
if m ~= n
disp('Matrix is must be square!')
end
nul = zeros(n);
bj = eye(n); bj1 = eye(n);
j = 0;
while logr|logs
j = j+1;
c = a * bj;
aj = -simplify(trace(c)/j);
bj1 = bj;
bj = simplify(a * bj1 + aj * eye(n));
if aj ~= 0
s = j; as = aj; bs1 = bj1;
logs = 0;
end
if bj == nul
r = j; logr = 0;
end
end
k = r-s;
temp = -as^(-k-1) * a^k * bs1^(k+1);
DrazinInverse2a = simplify(temp);
DrazinInverse2a = simplify(DrazinInverse2a);