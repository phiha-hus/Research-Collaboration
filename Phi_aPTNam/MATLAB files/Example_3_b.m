% Try linear programing
Example_3_v2

n = 3;
I = eye(n);

A_neq = [kron((P'-I),I) -reshape(I,n^2,1)]

b_neq = reshape(Abar,n^2,1)


for i = 1:100

f = rand(n^2+1,1)
%f = ones(n^2+1,1)
%f = zeros(n^2+1,1)

x = linprog(f,A_neq,b_neq)

D = reshape(x(1:end-1),n,n)

H = Abar + D * (I-P)

eig(H)

if max(real(eig(H)))<0 && norm(eig(H))<1e+2
    Abar - H * P
    break
end

end