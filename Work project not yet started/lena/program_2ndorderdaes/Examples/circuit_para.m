%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% initial values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t0,y0,y1] = circuit_para(type)

I0 = 2;
V0 = 1;
R = 2;
L = 1/3;
C = 1.5;

t0=0;
     
lambda = -R/(2*L)+sqrt(R^2/((2*L)^2)-1/(L*C));


y0=[L*C*lambda^2*exp(lambda*t0)
    exp(lambda*t0)
    R*C*lambda*exp(lambda*t0)  ];

if(strcmp(type,'bdf'))
 y1=[L*C*lambda^3*exp(lambda*t0)
    lambda*exp(lambda*t0)
    R*C*lambda^2*exp(lambda*t0)];
else
 y1=[L*C*lambda^3*exp(lambda*t0)
    lambda*exp(lambda*t0)
    R*C*lambda^2*exp(lambda*t0)];
end


    
    