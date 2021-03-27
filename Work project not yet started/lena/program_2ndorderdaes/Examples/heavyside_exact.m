function ye = heavyside_exact(t)

n=3;

f1=t^n*2*heaviside(t)-1;
df1=2*n*t^(n-1)*heaviside(t)+2*t^n*dirac(t);
ddf1=2*n*(n-1)*t^(n-2)*heaviside(t)+4*n*t^(n-1)*dirac(t);
dddf1=2*n*(n-1)*t^(n-2)*heaviside(t)+4*n*t^(n-1)*dirac(t);   % ffffffffff

f2=exp(t);
df2=exp(t);
ddf2=exp(t);
dddf2=exp(t);


ye=[f1+2*df1+t*ddf1-ddf2
    f2-(t+1)*f1-2*df1-t*ddf1+ddf2
    df1+2*ddf1+t*dddf1+ddf1-dddf2
    df2-(t+1)*df1-f1-2*ddf1-t*dddf1-ddf1+dddf2];


