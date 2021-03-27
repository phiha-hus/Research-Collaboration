
function ye = circuit_exact(t)

%ye=[-0.7*exp(-0.354*t)+2.7*exp(-5.646*t)];
R = 2;
L = 1/3;
C = 1.5;
lambda = -R/(2*L)+sqrt(R^2/((2*L)^2)-1/(L*C));

ye=[L*C*lambda^2*exp(lambda*t)
    exp(lambda*t)
    R*C*lambda*exp(lambda*t)  ];