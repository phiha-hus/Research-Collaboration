%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% initial values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t0,y0,y1] = tanhrot_para(type)

t0=-1;
n=3;
g=1;


f1=t0^n*tanh(g*t0+1)-1;
df1=n*t0^(n-1)*tanh(g*t0+1)+t0^n*(g-g*tanh(g*t0+1)^2);
ddf1=n*(n-1)*t0^(n-2)*tanh(g*t0+1)+2*n*t0^(n-1)*(g-g*tanh(g*t0+1)^2)-2*g*t0^n*tanh(g*t0+1)*(g-g*tanh(g*t0+1)^2);
dddf1=n*(n-1)*(n-2)*t0^(n-3)*tanh(g*t0+1)+3*n*(n-1)*t0^(n-2)*(g-g*tanh(g*t0+1)^2)-6*n*g*t0^(n-1)*tanh(g*t0+1)*(g-g*tanh(g*t0+1)^2)-2*g*t0^n*(g-g*tanh(g*t0+1)^2)^2+4*g^2*t0^n*tanh(g*t0+1)^2*(g-g*tanh(g*t0+1)^2);



f2=exp(t0);
df2=exp(t0);
ddf2=exp(t0);
dddf2=exp(t0);


y0=[f1+2*df1+t0*ddf1-ddf2
    f2-(t0+1)*f1-2*df1-t0*ddf1+ddf2];

     
if(strcmp(type,'bdf'))
  y1=[df1+2*ddf1+t0*dddf1+ddf1-dddf2
      df2-(t0+1)*df1-f1-2*ddf1-t0*dddf1-ddf1+dddf2];
else
  y1=[df1+2*ddf1+t0*dddf1+ddf1-dddf2
      df2-(t0+1)*df1-f1-2*ddf1-t0*dddf1-ddf1+dddf2];
end