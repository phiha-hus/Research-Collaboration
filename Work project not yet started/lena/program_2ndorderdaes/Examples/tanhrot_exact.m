function ye = tanhrot_exact(t)

n=3;
g=1;

f1=t^n*tanh(g*t+1)-1;
df1=n*t^(n-1)*tanh(g*t+1)+t^n*(g-g*tanh(g*t+1)^2);
ddf1=n*(n-1)*t^(n-2)*tanh(g*t+1)+2*n*t^(n-1)*(g-g*tanh(g*t+1)^2)-2*g*t^n*tanh(g*t+1)*(g-g*tanh(g*t+1)^2);
dddf1=n*(n-1)*(n-2)*t^(n-3)*tanh(g*t+1)+3*n*(n-1)*t^(n-2)*(g-g*tanh(g*t+1)^2)-6*n*g*t^(n-1)*tanh(g*t+1)*(g-g*tanh(g*t+1)^2)-2*g*t^n*(g-g*tanh(g*t+1)^2)^2+4*g^2*t^n*tanh(g*t+1)^2*(g-g*tanh(g*t+1)^2);

f2=exp(t);
df2=exp(t);
ddf2=exp(t);
dddf2=exp(t);


ye=[f1+2*df1+t*ddf1-ddf2
    f2-(t+1)*f1-2*df1-t*ddf1+ddf2
    df1+2*ddf1+t*dddf1+ddf1-dddf2
    df2-(t+1)*df1-f1-2*ddf1-t*dddf1-ddf1+dddf2];
