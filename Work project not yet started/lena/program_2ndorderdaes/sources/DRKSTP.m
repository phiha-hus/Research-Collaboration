
function [T,H,HOLD,IDID]=DRKSTP(EDIF, ADIF, FDIF, N, N3, T, H, HOLD, HMIN,ERRACC, UROUND, SAFE, FACL, FACR, QUOT1, QUOT2,M, ID, IA, IU, IREQ, LMAX, NSING, JSTART, LAST, X, XPRIME, E0, LDE0, A0, LDA0, E, LDE, A, LDA,EQ, LDEQ, AQ, LDAQ, Z, LDZ, ZQ, LDZQ, AH, LDAH, W, LDW, CONT, LDCONT, WT, Z1, Z2, Z3,Y, B, F, FQ, IPAR, RPAR, IWORK, WORK, LWORK, PRED, MCONST)

%     Performs one step  of the Runge-Kutta integration of GELDA [2].
%     ******************************************************************

ZERO = 0.0;
ONE = 1.0;

ONST   = 11;
ONEV   = 12;
ONFA   = 13;
OETF   = 14;
OCTF   = 15;
OIPVT  = 21;

% *** *** *** *** *** *** ***
%  INITIALISATIONS
% *** *** *** *** *** *** ***
if (JSTART == 0)
    IWORK(OETF) = 0;
    IWORK(OCTF) = 0;
    NSING = 0;
end

% ---------- Constants ---------
      SQ6=sqrt(6.0);
      C1=(4.D0-SQ6)/10.D0;
      C2=(4.D0+SQ6)/10.D0;
      C1M1=C1-1.D0;
      C2M1=C2-1.D0;
      C1MC2=C1-C2;
      DD1=-(13.D0+7.D0*SQ6)/3.D0;
      DD2=(-13.D0+7.D0*SQ6)/3.D0;
      DD3=-1.D0/3.D0;
      U1=(6.0+81.0^(1.0/3.0)-9.0^(1.D0/3.D0))/30.D0;
      U1=ONE/U1;
      A11=( 88.D0-  7.D0*SQ6)/360.D0;
      A12=(296.D0-169.D0*SQ6)/ 18.D2;
      A13=( -2.D0+  3.D0*SQ6)/225.D0;
      A21=(296.D0+169.D0*SQ6)/ 18.D2;
      A22=( 88.D0+  7.D0*SQ6)/360.D0;
      A23=( -2.D0-  3.D0*SQ6)/225.D0;
      A31=( 16.D0-       SQ6)/ 36.D0;
      A32=( 16.D0+       SQ6)/ 36.D0;
      A33= 1.D0/9.D0;
%***
IDID=0;
REJECT=false;
if (JSTART>0 && MCONST && IU==0 && H==HOLD)
    FACTOR=false;
else
    FACTOR=true;
end
N2=2*N;
HHFAC=H;
if (~MCONST || JSTART == 0)
    for  i=1:N
        for j=1:N
            E0(i,j)=E(i,j);
            A0(i,j)=A(i,j);
        end
    end
end



% *** *** *** *** *** *** ***
%  BASIC INTEGRATION STEP
% *** *** *** *** *** *** ***
%  10  CONTINUE
CONTINUE=true;

while(CONTINUE)
    if (abs(H)<=HMIN)
        IDID=-6;
        return;
    end
    TPH=T+H;
    HA11=H*A11;
    HA12=H*A12;
    HA13=H*A13;
    HA21=H*A21;
    HA22=H*A22;
    HA23=H*A23;
    HA31=H*A31;
    HA32=H*A32;
    HA33=H*A33;

    % --- Evaluate E,A and F with respect to the three RK stages
    % ---   and compute solution matrix and right-hand side

    % --- First Runge-Kutta stage
    DNFIX1(EDIF, ADIF, FDIF, N, T+C1*H, LMAX, M, ID, IA, IU, IREQ, MQ, IDQ, IAQ, IUQ, E, LDE, A, LDA, F, EQ, LDEQ, AQ, LDAQ, FQ, Z, LDZ, ZQ, LDZQ,AH, LDAH, IPAR, RPAR, WORK, LWORK, MCONST, IER)
    if(IER ~= 0)
        IDID=-10;
        return;
    end
    IWORK(ONEV)=IWORK(ONEV)+1;
    DGEMV('N',N,N,ONE,A,LDA,X,1,ONE,F,1);
    DCOPY(N,F,1,B(1),1);
    if ((~MCONST) && FACTOR)
        for i=1:N
            for j=1:N
                W(i,j   )=E(i,j) -HA11*A(i,j);
                W(i,N+j )=       -HA12*A(i,j);
                W(i,N2+j)=       -HA13*A(i,j);
            end
        end
    end

    % --- Second Runge-Kutta stage
    DNFIX1(EDIF, ADIF, FDIF, N, T+C2*H, LMAX, M, ID, IA, IU, IREQ, MQ, IDQ, IAQ, IUQ, E, LDE, A, LDA, F, EQ, LDEQ, AQ, LDAQ, FQ, Z, LDZ, ZQ, LDZQ,AH, LDAH, IPAR, RPAR, WORK, LWORK, MCONST, IER)
    if (IER ~= 0)
        IDID=-10;
        return;
    end
    IWORK(ONEV)=IWORK(ONEV)+1;
    DGEMV('N',N,N,ONE,A,LDA,X,1,ONE,F,1);
    DCOPY(N,F,1,B(N+1),1);
    if ((~ MCONST) && FACTOR)
        for i=1:N
            for j=1:N
                W(N+i,j   )=       -HA21*A(i,j);
                W(N+i,N+j )=E(i,j) -HA22*A(i,j);
                W(N+i,N2+j)=       -HA23*A(i,j);
            end
        end
    end

    % --- Third Runge-Kutta stage
    DNFIX1(EDIF, ADIF, FDIF, N, TPH, LMAX, M, ID, IA, IU,IREQ, MQ, IDQ, IAQ, IUQ, E, LDE, A, LDA, F,EQ, LDEQ, AQ, LDAQ, FQ, Z, LDZ, ZQ, LDZQ,AH, LDAH, IPAR, RPAR, WORK, LWORK, MCONST, IER)
    if (IER ~= 0)
        IDID=-10;
        return;
    end
    IWORK(ONEV)=IWORK(ONEV)+1;
    DGEMV('N',N,N,ONE,A,LDA,X,1,ONE,F,1);
    DCOPY(N,F,1,B(N2+1),1);
    if ((~ MCONST) && FACTOR)
        for i=1:N
            for j=1:N
                W(N2+i,j   )=       -HA31*A(i,j);
                W(N2+i,N+j )=       -HA32*A(i,j);
                W(N2+i,N2+j)=E(i,j) -HA33*A(i,j);
            end
        end
    end

    % --- Decomposition and solution of the Runge-Kutta system
    FAC1=U1/H;
    if (MCONST)
        DRKCXS(N,N2,N3,IU,FAC1,H,UROUND,AH,LDAH,W,N, E0,LDE0,A0,LDA0,B,Y,IWORK,WORK,FACTOR,IER);
    else
        DRKRLS(N,N3,IU,FAC1,UROUND,AH,LDAH,W,LDW, E0,LDE0,A0,LDA0,B,Y,IWORK,WORK,IER);
    end
    if (IER~=0)
        if (IER~=0)
            NSING=NSING+1;
            if(NSING>=5)
                IDID=-8;
                return;
            end
        end
        H=H*0.50;
        HHFAC=0.50;
        REJECT=true;
        LAST=false;
        FACTOR=true;
        continue;
    end

    % --- Compute the vectors Z1, Z2, Z3
    DSCAL(N,ZERO,Z1,1);
    DAXPY(N,HA11,Y      ,1,Z1,1);
    DAXPY(N,HA12,Y(N+1) ,1,Z1,1);
    DAXPY(N,HA13,Y(N2+1),1,Z1,1);
    DSCAL(N,ZERO,Z2,1);
    DAXPY(N,HA21,Y      ,1,Z2,1);
    DAXPY(N,HA22,Y(N+1) ,1,Z2,1);
    DAXPY(N,HA23,Y(N2+1),1,Z2,1);
    DSCAL(N,ZERO,Z3,1);
    DAXPY(N,HA31,Y      ,1,Z3,1);
    DAXPY(N,HA32,Y(N+1) ,1,Z3,1);
    DAXPY(N,HA33,Y(N2+1),1,Z3,1);

    % *** *** *** *** *** *** ***
    % ERROR ESTIMATION
    % *** *** *** *** *** *** ***
    HEE1=DD1/H;
    HEE2=DD2/H;
    HEE3=DD3/H;
    DCOPY(N,XPRIME,1,F,1);
    DAXPY(N,HEE1,Z1,1,F,1);
    DAXPY(N,HEE2,Z2,1,F,1);
    DAXPY(N,HEE3,Z3,1,F,1);

    % --- Consider enlarged system
    DGEMV('N',N,N,ONE,E0,LDE0,F,1,ZERO,CONT(1,1),1);
    DGEMV('N',N,N,ONE,A0,LDA0,F,1,ZERO,CONT(1,2),1);
    if (IU == 0)
        DGETRS ('No transpose',N,2,AH,LDAH,IWORK(OIPVT+N3),CONT,LDCONT,IER);
        if (IER ~= 0)
            IDID=-6;
            return;
        end
    else
        for i=1:N
            for j=1:N
                AH(i,j) = FAC1*E0(i,j) - A0(i,j);
            end
            IWORK(OIPVT+N3+i-1)=0;
        end
        DGELSX(N,N,2,AH,LDAH,CONT,LDCONT,IWORK(OIPVT+N3),UROUND,IRANK,WORK,IER);
        if(IER ~= 0)
            IDID=-6;
            return;
        end
    end

    % --- The error norm (Index-2-components must be multiplied by H)
    ERR=ZERO;
    for i=1:N
        ERR=ERR + (CONT(i,1)*WT(i))^2 + (CONT(i,2)*HHFAC*WT(i))^2;
    end
    ERR=MAX(SQRT(ERR/N2),1.D-10);
    if (ERR >= ONE && (JSTART == 0 || REJECT))

        % ---    Use improved error estimator
        DAXPY(N,ONE,CONT(1,2),1,F,1);
        DGEMV('N',N,N,ONE,E0,LDE0,F,1,ZERO,CONT(1,1),1);
        DGEMV('N',N,N,ONE,A0,LDA0,F,1,ZERO,CONT(1,2),1);
        if (IU == 0)
            DGETRS ('No transpose',N,2,AH,LDAH,IWORK(OIPVT+N3),CONT,LDCONT,IER);
            if (IER ~= 0)
                IDID=-6;
                return;
            end
        else
            for i=1:N
                for j=1:N
                    AH(i,j) = FAC1*E0(i,j) - A0(i,j);
                end
                IWORK(OIPVT+N3+i-1)=0;
            end
            DGELSX(N,N,2,AH,LDAH,CONT,N,IWORK(OIPVT+N3),UROUND,IRANK,WORK,IER);
            if (IER ~= 0)
                IDID=-6;
                return;
            end
        end

        % ---    The error norm (Index-2-components must be multiplied by H)
        ERR=ZERO;
        for i=1:N
            ERR=ERR + (CONT(i,1)*WT(i))^2 + (CONT(i,2)*HHFAC*WT(I))^2;
        end
        ERR=MAX(SQRT(ERR/N2),1.D-10);
    end

    % --- Computation of HNEW (we require 0.2 <= HNEW/H <= 8.0)
    FAC=SAFE;
    QUOT=MAX(FACR,MIN(FACL,ERR^.25D0/FAC));
    HNEW=H/QUOT;

    % --- Is the error small enough ?
    if (ERR<ONE)
        % *** *** *** *** *** *** ***
        % STEP IS ACCEPTED
        % *** *** *** *** *** *** ***
        if (JSTART == 0)
            JSTART=1;
        end
        IWORK(ONST)=IWORK(ONST)+1;
        if (PRED)
            % ---    Predictive controller of Gustafsson
            if (IWORK(ONST)>1)
                FACGUS=(HOLD/H)*(ERR^2/ERRACC)^0.25D0/SAFE;
                FACGUS=MAX(FACR,MIN(FACL,FACGUS));
                QUOT=MAX(QUOT,FACGUS);
                HNEW=H/QUOT;
            end
            ERRACC=MAX(1.0D-2,ERR);
        end
        HOLD=H;

        % ---    Update T, solution vector X and derivative vector XPRIME
        T=TPH;
        DAXPY(N,ONE,Z3,1,X,1);
        DCOPY(N,Y(N2+1),1,XPRIME,1);
        if(LAST)
            return;
        end
        if (REJECT)
            POSNEG=SIGN(ONE,H);
            HNEW=POSNEG*MIN(ABS(HNEW),ABS(H));
        end
        QT=HNEW/H;
        if (~MCONST || IU~=0 || QT<QUOT1 ||  QT>QUOT2)
            H=HNEW;
        end
        return;
    else
        % *** *** *** *** *** *** ***
        % STEP IS REJECTED
        % *** *** *** *** *** *** ***
        REJECT=true;
        LAST=false;
        FACTOR=true;
        if (JSTART == 0)
            H=H*0.10;
            HHFAC=0.10;
        else
            HHFAC=HNEW/H;
            H=HNEW;
        end
        IWORK(OETF)=IWORK(OETF)+1;
        continue;
    end
end

