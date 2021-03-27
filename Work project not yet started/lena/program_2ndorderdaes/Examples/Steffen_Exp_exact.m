function ye = Steffen_Exp_exact(t)

a = 10;

ye = [exp(-a*t);-exp(-a*t);2*exp(-a*t)];