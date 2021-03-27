function ye = circuit_exact_1_start(t0,h)

%ye=[-0.7*exp(-0.354*t)+2.7*exp(-5.646*t)];
R = 2;
L = 1/3;
C = 1.5;
lambda = -R/(2*L)+sqrt(R^2/((2*L)^2)-1/(L*C));

ye=[L*C*lambda^2*exp(lambda*t0)
    exp(lambda*t0)
    R*C*lambda*exp(lambda*t0) 
    C*lambda*exp(lambda*t0)
    h*L*C*lambda^3*exp(lambda*t0)
    h*lambda*exp(lambda*t0)
    h*R*C*lambda^2*exp(lambda*t0) 
    h*C*lambda^2*exp(lambda*t0)
    h^2*L*C*lambda^4*exp(lambda*t0)
    h^2*lambda^2*exp(lambda*t0)
    h^2*R*C*lambda^3*exp(lambda*t0) 
    h^2*C*lambda^3*exp(lambda*t0)];