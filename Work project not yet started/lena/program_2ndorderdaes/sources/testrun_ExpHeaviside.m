% testrun_ExpHeaviside

%     .. Data Statements ..
DTOUT =[0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00];

NEQ =2; 
LRW = 5000; 
LIW = 100;
LUN = 6;
NOUT = 10;
TSTART = 0.0;

METHOD = 1; % BDF

for i=1:20
         INFO(i) = 0;
end

INFO(8)=1;
RWORK(3)=0.001;

% starting time and consistent initial values
T = TSTART;
X(1)=0.0;
X(2)=0.0;
XPRIME(1)=0.0D0;
XPRIME(2)=0.0D0;

% tolerances
ATOL(1) = 1.0e-5;
RTOL(1) = 1.0e-5;

F='heaviside';
dF='heaviside_d1';
ddF='heaviside_d2';


% solve the problem.
time = clock;

for IOUT = 1:NOUT
         TOUT = DTOUT(IOUT);
         [T,X,XPRIME,H,RWORK,IWORK,IERR]=HOBDF(F,dF,ddF,NEQ,T,TOUT,X,XPRIME,RTOL,ATOL,METHOD,INFO,LRW,LIW,RWORK);
         if(IERR < 0)
             error(' error!!!');
         end
         %ERROR(1) = ((1+T)*exp(-T) - X(1)) /(abs((1+T)*exp(-T))+1);
         %ERROR(2) = (exp(-T) - X(2)) /(exp(-T)+1);
         %ERO = max(ERO,norm(NEQ,ERROR,1)/sqrt(NEQ));
         %HU = RWORK(7);
         %NQU = IWORK(8);
end

%     Finally, we display some final statistics
etime(clock, time)
NST = IWORK(11);
NFE = IWORK(12);
NJE = IWORK(13);
