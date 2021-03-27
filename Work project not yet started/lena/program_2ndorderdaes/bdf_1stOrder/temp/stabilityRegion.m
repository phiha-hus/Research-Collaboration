% stability matrixx
function  stabilityRegion(A,At,U,Ut,Vt,Bt,n)

c =1;

for x = -5:0.01:0
    for y = -5:0.01:5
        z = x+i*y;
        M = Vt+z*Bt*inv(eye(n,n)-z*At)*Ut;
        if( norm(M,inf)<=1 )
            list(c)=z;
            c = c+1;
            %eig(M)
        end
    end
end

for k=1:c-1
    plot(list(k));
    hold on;
end