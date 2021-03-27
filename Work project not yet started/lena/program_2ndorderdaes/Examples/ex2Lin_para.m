%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% initial values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t0,y0,y1] = ex2Lin_para(type)

t0=-1;


f1=sin(t0)+t0^2*cos(t0);
df1=(2*t0+1)*cos(t0)-t0^2*sin(t0);
ddf1=(2-t0^2)*cos(t0)-(4*t0+1)*sin(t0);
dddf1=(-6*t0-1)*cos(t0)-(6-t0^2)*sin(t0);  

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