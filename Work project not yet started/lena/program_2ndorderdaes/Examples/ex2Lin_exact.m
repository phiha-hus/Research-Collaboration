function ye = ex2Lin_exact(t)

f1=sin(t)+t^2*cos(t);
df1=(2*t+1)*cos(t)-t^2*sin(t);
ddf1=(2-t^2)*cos(t)-(4*t+1)*sin(t);
dddf1=(-6*t-1)*cos(t)-(6-t^2)*sin(t);  

f2=exp(t);
df2=exp(t);
ddf2=exp(t);
dddf2=exp(t);


ye=[f1+2*df1+t*ddf1-ddf2
    f2-(t+1)*f1-2*df1-t*ddf1+ddf2
    df1+2*ddf1+t*dddf1+ddf1-dddf2
    df2-(t+1)*df1-f1-2*ddf1-t*dddf1-ddf1+dddf2];