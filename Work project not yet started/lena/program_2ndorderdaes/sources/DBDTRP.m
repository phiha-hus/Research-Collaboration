
function DBDTRP(T, TOUT, NEQ, KOLD, PHI, PSI, XOUT, XPOUT)

%     Interpolation routine for DBDFST.

%
%     METHOD
%
%     The methods in subroutine DBDSTP use polynomials to approximate
%     the solution. DBDTRP approximates the solution and its
%     derivative at time TOUT by evaluating one of these polynomials,
%     and its derivative, there. Information defining this polynomial
%     is passed from DBDSTP, so DBDTRP cannot be used alone.
%
%     DBDTRP is a modified version of the DDATRP subroutine of DASSL [3].

%  ******************************************************************


ZERO = 0.0; 
ONE = 1.0;

KOLDP1=KOLD+1;
TEMP1=TOUT-T;

DCOPY(NEQ,PHI,1,XOUT,1);
for j=1:NEQ
    XPOUT(j) = ZERO;
end
C=ONE;
D=ZERO;
GAMMA=TEMP1/PSI(1);
for j=2:KOLDP1
    D=D*GAMMA+C/PSI(j-1);
    C=C*GAMMA;
    GAMMA=(TEMP1+PSI(j-1))/PSI(j);

    DAXPY(NEQ,C,PHI(1,j),1,XOUT,1);
    DAXPY(NEQ,D,PHI(1,j),1,XPOUT,1);
end
