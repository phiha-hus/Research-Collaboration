% output quantities
function    [y] = glm_external( yn,Ypp,n,M,h,i,s,k )

Bt = M(2*s+1:2*s+s,1:s );
Vt = M(2*s+1:2*s+s,s+1:s+k );

if( i>1)
    D = zeros(k,k);
    for j=1:k
        D(j,j)=(h(i)/h(i-1))^(j-1);
    end
    
    yn = kron(D,eye(n,n))*yn;
end

y = kron(Vt,eye(n,n))*yn + h(i)^2*kron(Bt,eye(n,n))*Ypp;
