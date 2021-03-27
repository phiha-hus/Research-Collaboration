% output quantities
function    [y] = glm1_external( yn,Yp,n,M,h,i,p,s,k )

B = M(s+1:s+s,1:s );
V = M(s+1:s+s,s+1:s+k );

if( i>1)
    D = zeros(p+1,p+1);
    for j=1:p+1
        D(j,j)=(h(i)/h(i-1))^(j-1);
    end
    
    yn = kron(D,eye(n,n))*yn;
end

y = kron(V,eye(n,n))*yn + h(i)*kron(B,eye(n,n))*Yp;
