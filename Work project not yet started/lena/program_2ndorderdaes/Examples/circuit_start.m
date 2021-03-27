function ye = circuit_start(t0,h,s)

R = 2;
L = 1/3;
C = 1.5;
lambda = -R/(2*L)+sqrt(R^2/((2*L)^2)-1/(L*C));


if( s==3)
ye=[L*C*lambda^2*exp(lambda*t0)
    exp(lambda*t0)
    R*C*lambda*exp(lambda*t0) 
    h*L*C*lambda^3*exp(lambda*t0)
    h*lambda*exp(lambda*t0)
    h*R*C*lambda^2*exp(lambda*t0)  
    h^2*L*C*lambda^4*exp(lambda*t0)
    h^2*lambda^2*exp(lambda*t0)
    h^2*R*C*lambda^3*exp(lambda*t0) ];

end

if( s==4)

ye=[L*C*lambda^2*exp(lambda*t0)
    exp(lambda*t0)
    R*C*lambda*exp(lambda*t0) 
    h*L*C*lambda^3*exp(lambda*t0)
    h*lambda*exp(lambda*t0)
    h*R*C*lambda^2*exp(lambda*t0)  
    h^2*L*C*lambda^4*exp(lambda*t0)
    h^2*lambda^2*exp(lambda*t0)
    h^2*R*C*lambda^3*exp(lambda*t0) 
    h^3*L*C*lambda^5*exp(lambda*t0)
    h^3*lambda^3*exp(lambda*t0)
    h^3*R*C*lambda^4*exp(lambda*t0)];

end

if( s==2)
ye=[L*C*lambda^2*exp(lambda*t0)
    exp(lambda*t0)
    R*C*lambda*exp(lambda*t0) 
    h*L*C*lambda^3*exp(lambda*t0)
    h*lambda*exp(lambda*t0)
    h*R*C*lambda^2*exp(lambda*t0)   ];

end

if( s==1)
ye=[L*C*lambda^2*exp(lambda*t0)
    exp(lambda*t0)
    R*C*lambda*exp(lambda*t0)    ];

end