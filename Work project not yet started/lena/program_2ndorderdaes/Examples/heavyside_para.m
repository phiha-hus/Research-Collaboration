%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% initial values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t0,y0,y1] = heavyside_para(type)

t0=-1;
n=3;

f1=t0^n*2*heaviside(t0)-1;
df1=2*n*t0^(n-1)*heaviside(t0)+2*t0^n*dirac(t0);
ddf1=2*n*(n-1)*t0^(n-2)*heaviside(t0)+4*n*t0^(n-1)*dirac(t0);
dddf1=2*n*(n-1)*(n-2)*t0^(n-3)*heaviside(t0)+6*n*(n-1)*t0^(n-2)*dirac(t0);

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