function ye = Steffen_Exp_start(t0,h)

a = 10;

ye = [exp(-a*t0);
     -exp(-a*t0);
     2*exp(-a*t0)
     -h*a*exp(-a*t0);
     h*a*exp(-a*t0);
     -h*2*a*exp(-a*t0)
     h^2*a^2*exp(-a*t0);
     -h^2*a^2*exp(-a*t0);
     h^2*2*a^2*exp(-a*t0)];