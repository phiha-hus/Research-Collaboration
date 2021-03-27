function []=BDFSTP (FDIF, NEQ, T, H, HOLD, HMIN, K, KOLD, M, IREQ, LMAX,NS, IPHASE, ROWEQU, COLEQU, JSTART, X, XPRIME, E, LDE, A, LDA, Z1, LDZ1, Z2Q, LDZ2Q, WT, W,LDW, PHI, LDPHI, ALPHA, BETA, GAMMA, PSI, SIGMA, DELTA, ERRV, RS, CS, F, EQ, LDEQ, AQ, LDAQ, FQ, AH, LDAH, IPAR, RPAR, IWM, WORK, LWORK, MCONST, IWARN, IDID)

%     Performs one step of the BDF integration of GELDA [2].
%     DBDSTP solves a system of differential-algebraic equations of
%     the form E(T) XPRIME(T) = A(T) X(T) + F(T), for one step
%     (normally from T to T+H).

%{  
function [ y,converged, IER ] = bdf2_step(funct,funct_d1,funct_d2,funct_d3,n,nL,KK,TT,t,h,i,y,RTOL,ATOL,exponent )

    % BDF coefficients of BDF(l,j,k)
    a = bdfk_coef(KK,TT);
    
    nq = n-nL;  % nL=1 !!

    %predictor value
    sum1=0;
    for ii=1:KK(1)+1
        sum1=sum1+phi_star(ii,i,y,t(1:i),h);
    end

    %starting value for Newton-Iteration
    x(:,1)=sum1;
    d(:,1)=sum1;

    % Newton-Iteration
    m=1;          
    maxIt = 4; 
    abbruch=1;
    converged = 0;

    while(abbruch >= 0.33 && m<=maxIt)
    
        y_0=y(1:n,:);
        y_1=y(n+1:n+nq,:);
    
        % Zusammenbau von F
        sum1 = 0;
        sum2 = 0;
        sum3 = 0;
        term_a = 0;
    
        for c=1:KK(1)
            sum1 = sum1 + a(c)*divDiff( t(i+1-c:i+1),[y_1(:,i+1-c:i),x(n+1:n+nq,m)],c);
            prod=1;
            for j=1:c-1
                prod = prod*(t(i+1)-t(i+1-j));
            end
            sum3 = sum3 + prod*divDiff(t(i+1-c:i+1),[y_0(1:nq,i+1-c:i),x(1:nq,m)],c);
            sum2 = sum2 + 1/(t(i+1)-t(i+1-c));
            prod2=1;
            for j=1:c
                prod2=prod2*(t(i+1)-t(i+1-j));
            end
            term_a=term_a+a(c)*1/prod2;
        end
        
        F =  [feval(funct,    t(i+1),x(1:n,m),x(n+1:n+nq,m), sum1 )
              x(n+1:n+nq,m)-sum3];
              
              
        D1=feval(funct_d1,t(i+1),x(1:n,m),x(n+1:n+nq,m), sum1);
        D2=feval(funct_d3,t(i+1),x(1:n,m),x(n+1:n+nq,m), sum1,'bdf')*term_a-feval(funct_d2,t(i+1),x(1:n,m),x(n+1:n+nq,m), sum1,'bdf');
        D3= eye(nq,nq+nL).*(-sum2);
        D4= eye(nq,nq);
        DF= [D1 D2
             D3 D4];

        if( det(DF) == 0) % J singular
            IER=-1;
        else
            IER=0; % J regular
        end; 
    
        % Gauss Zerlegung
        [L,U,P]=lu(DF);
        bb = P*(-F);
        %yy = L\bb;
        %d(:,m+1) = U\yy;
        yy = forward(L,bb);
        d(:,m+1) = backward(U,yy);
        
        
        x(:,m+1)=x(:,m)+d(:,m+1);
    
        m=m+1;
        if(m>2)
            rho = (myNorm(x(:,m)-x(:,m-1),RTOL,ATOL)/myNorm(x(:,2)-x(:,1),RTOL,ATOL))^(1/(m-1));
            abbruch = rho/(1-rho)*myNorm(x(:,m)-x(:,m-1),RTOL,ATOL);
            if(rho>0.9)
                converged = 1;
                y=x(:,m);
                return;
            end
        end
    end

    if(m > maxIt)
        converged = 1;
    end

    y=x(:,m);

%}

      ONE=1.0;
      OMXORD=3; 
      ONST=11; 
      ONEV=12; 
      ONFA=13; 
      OETF=14; 
      OCTF=15; 
      OIPVT=21;

%----------------------------------------------------------------------
%     Block 1.
%     Initialize. On the first call, set
%     the order to 1 and initialize
%     other variables.
%----------------------------------------------------------------------

%     Initializations for all calls.
      IDID=1;
      TOLD=T;
      NCF=0;
      NSF=0;
      NEF=0;
      IRED=0;

%     If this is the first step, perform  other initializations.
      if(JSTART ~= 1) 
        K=1;
        KOLD=0;
        HOLD=0.00;
        PSI(1)=H;
        %CJOLD = 0.00;
        FACTOR = true;
        IPHASE = 0;
        NS=0;
        JSTART=1;
      end

%----------------------------------------------------------------------
%     Block 2
%     Compute coefficients of formulas for
%     this step.
%----------------------------------------------------------------------

 %200  CONTINUE
      KP1=K+1;
      KP2=K+2;
      KM1=K-1;
      TOLD=T;
      if (H~=HOLD || K ~= KOLD) 
          NS = 0;
      end
      NS=MIN(NS+1,KOLD+2);
      NSP1=NS+1;

%{      
      if (KP1 >= NS) 
        BETA(1)=1.0D0;
        ALPHA(1)=1.0D0;
        TEMP1=H;
        GAMMA(1)=0.0D0;
        SIGMA(1)=1.0D0;
        for I=2:KP1
           TEMP2=PSI(I-1);
           PSI(I-1)=TEMP1;
           BETA(I)=BETA(I-1)*PSI(I-1)/TEMP2;
           TEMP1=TEMP2+H;
           ALPHA(I)=H/TEMP1;
           SIGMA(I)=(I-1)*SIGMA(I-1)*ALPHA(I);
           GAMMA(I)=GAMMA(I-1)+ALPHA(I-1)/H;
        end
        PSI(KP1)=TEMP1;
      end

%     Compute ALPHAS, ALPHA0.
      ALPHAS = 0.00;
      for I=1:K
        ALPHAS = ALPHAS - 1.0D0/I;
      end
      ALPHA0 = 0.0D0;
      for I = 1:K
        ALPHA0 = ALPHA0 - ALPHA(I);
      end

%     Compute leading coefficient CJ
%     NOTE: If MCONST = TRUE and IU = 0 we can use the old factorization
      CJ = -ALPHAS/H;
      if (MCONST && IU==0) 
        if (CJ==CJOLD) 
          FACTOR = false;
        else
          FACTOR = true;
        end
      end

%     Compute variable stepsize error coefficient CK.
      CK = ABS(ALPHA(KP1) + ALPHAS - ALPHA0);
      CK = MAX(CK,ALPHA(KP1));

%     Change PHI to PHI STAR.
      if (KP1 >= NSP1)
        for J=NSP1:KP1
            DSCAL(NEQ,BETA(J),PHI(1,J),1)
        end
end
%}

%     Update time.
      T=T+H;
      
%{
%----------------------------------------------------------------------
%     Block 3
%     Predict the solution and derivative,
%     and solve the corrector equation.
%----------------------------------------------------------------------

%     First, predict the solution and derivative.
% 300  CONTINUE

      %DCOPY(NEQ,PHI(1,1),1,X,1)
      X=PHI(1,1);
      for I=1:NEQ
        XPRIME(I)=0.0D0
      end

      for J=2:KP1
        DAXPY(NEQ,ONE,PHI(1,J),1,X,1);
        DAXPY(NEQ,GAMMA(J),PHI(1,J),1,XPRIME,1);
      end

%     Compute E(T), A(T), F(T) and set up linear system.
      IWM(ONEV) = IWM(ONEV) + 1;
      DNFIX1(EDIF, ADIF, FDIF, NEQ, T, LMAX, M, ID, IA, IU, IREQ,  MQ, IDQ, IAQ, IUQ, E, LDE, A, LDA, F, EQ, LDEQ, AQ, LDAQ, FQ, Z1, LDZ1, Z2Q, LDZ2Q, AH, LDAH, IPAR, RPAR, WORK, LWORK, MCONST, IRED)
      if (IRED < 0) 
          GOTO 380
      end
      if (FACTOR) 
        for I = 1:NEQ
          for  J = 1:NEQ
            W(I,J) = -CJ*E(I,J) + A(I,J)
          end
        end
      end
      %DCOPY(NEQ,XPRIME,1,ERRV,1);
      ERRV=XPRIME;
      DAXPY(NEQ,-CJ,X,1,ERRV,1);
      %DCOPY(NEQ,F,1,DELTA,1);
      DELTA=F;
      DGEMV('N',NEQ,NEQ,1.0D0,E,LDE,ERRV,1,-1.0D0,DELTA,1);
C
C     First case: IU = 0 and we have to solve a linear system.
    
         IF (FACTOR) THEN
            CALL DGEEQU(NEQ,NEQ,W,NEQ,RS,CS,ROWCND,COLCND,AMAX,IER)
            CALL DLAQGE(NEQ,NEQ,W,NEQ,RS,CS,ROWCND,COLCND,AMAX,EQUED)
            ROW = (LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' ))
            IF (ROW) THEN
               ROWEQU = 1
            ELSE
               ROWEQU = 0
            ENDIF
            COL = (LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' ))
            IF (COL) THEN
               COLEQU = 1
            ELSE
               COLEQU = 0
            ENDIF
            CALL DGETRF(NEQ,NEQ,W,NEQ,IWM(OIPVT),IER)
            IF (IER .NE. 0) GOTO 380
            CJOLD = CJ
            IWM(ONFA) = IWM(ONFA) + 1
         ENDIF
         IF(ROWEQU.EQ.1) THEN
            DO 930 I = 1, NEQ
               DELTA(I) = RS(I)*DELTA(I)
 930        CONTINUE
         ENDIF
         CALL DGETRS('N',NEQ,1,W,NEQ,IWM(OIPVT),DELTA,NEQ,IER)
         IF(COLEQU.EQ.1) THEN
            DO 970 I = 1, NEQ
               DELTA(I) = CS(I)*DELTA(I)
 970        CONTINUE
         ENDIF
C
C           Update ERRV, X and XPRIME:
C           If solution of the system is stored in DELTA (not in ERRV)
C           we get.
C             ERRV   = DELTA - X
C             X      = DELTA
C             XPRIME = XPRIME * CJ * (DELTA - X)
         CALL DCOPY(NEQ,DELTA,1,ERRV,1)
         CALL DAXPY(NEQ,-1.0D0,X,1,ERRV,1)
         CALL DCOPY(NEQ,DELTA,1,X,1)
         CALL DAXPY(NEQ,CJ,ERRV,1,XPRIME,1)

 380  IF((IER .NE. 0) .OR. (IRED .LT. 0)) GOTO 600

C----------------------------------------------------------------------
C     Block 4
C     Estimate the errors at orders K,K-1,K-2
C     As if constant stepsize was used. Estimate
C     the local error at order K and test
C     whether the current step is successful.
C----------------------------------------------------------------------

C     Estimate errors at orders K,K-1,K-2.
      ENORM = DGENRM(NEQ,ERRV,WT)
      ERK = SIGMA(K+1)*ENORM
      TERK = (K+1)*ERK
      EST = ERK
      KNEW=K
      IF(K .NE. 1) THEN
C
C       Set DELTA(*) = PHI(*,KP1) + ERRV(*).
        DCOPY(NEQ,PHI(1,KP1),1,DELTA,1)
        DAXPY(NEQ,ONE,ERRV,1,DELTA,1)
        ERKM1=SIGMA(K)*DGENRM(NEQ,DELTA,WT)
        TERKM1 = K*ERKM1
        IF(K .LE. 2) THEN
          IF(TERKM1 .LE. 0.5D0*TERK) THEN
C
C           Lower the ordder.
            KNEW=K-1
            EST = ERKM1
	  ENDIF
        ELSE
C
C         Set DELTA(*) = PHI(*,K) + DELTA(*).
          CALL DAXPY(NEQ,ONE,PHI(1,K),1,DELTA,1)
          ERKM2=SIGMA(K-1)*DGENRM(NEQ,DELTA,WT)
          TERKM2 = (K-1)*ERKM2
          IF(MAX(TERKM1,TERKM2).LE.TERK) THEN
C
C           Lower the order.
            KNEW=K-1
            EST = ERKM1
	  ENDIF
	ENDIF
      ENDIF

%     Calculate the local error for the current step
%     to see if the step was successful.
      ERR = CK * ENORM
      if(ERR > 1.0D0)
        GO TO 600
      end

C----------------------------------------------------------------------
C     BLOCK 5
C     The step is successful. Determine the best order and stepsize for
C     the next step. Update the differences for the next step.
C----------------------------------------------------------------------
      IDID=1;
      IWM(ONST)=IWM(ONST)+1;
      KDIFF=K-KOLD;
      KOLD=K;
      HOLD=H;

%     Estimate the error at order K+1 unless:
%        Already decided to lower order, or already using maximum order, or
%        stepsize not constant, or order raised in previous step.
      IF (KNEW.EQ.KM1.OR.K.EQ.IWM(OMXORD)) IPHASE=1
      IF (IPHASE .NE. 0) THEN
        IF (KNEW.EQ.KM1) GOTO 540
        IF (K.EQ.IWM(OMXORD)) GOTO 550
        IF (KP1.GE.NS.OR.KDIFF.EQ.1) GOTO 550

C%       Set DELTA(*)=ERRV(*)-PHI(*,KP2).
        CALL DCOPY(NEQ,ERRV,1,DELTA,1)
        CALL DAXPY(NEQ,-ONE,PHI(1,KP2),1,DELTA,1)
        ERKP1 = (1.0D0/(K+2))*DGENRM(NEQ,DELTA,WT)
        TERKP1 = (K+2)*ERKP1
        IF(K.LE.1) THEN
          IF(TERKP1.GE.0.5D0*TERK)GO TO 550
	ELSE
          IF(TERKM1.LE.MIN(TERK,TERKP1))GO TO 540
          IF(TERKP1.GE.TERK.OR.K.EQ.IWM(OMXORD))GO TO 550
	ENDIF

%       Raise order.
        K=KP1
        EST = ERKP1
        GOTO 550

%       Lower order.
 540    K=KM1
        EST = ERKM1
        GOTO 550

%       Determine the appropriate stepsize for  the next step.
 550    HNEW=H
        TEMP2=K+1
        R=(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
        IF (R .GE. 2.0D0) THEN
          HNEW = 2.0D0*H
        ELSEIF(R .LE. 1.0D0) THEN
          R = MAX(0.5D0,MIN(0.9D0,R))
          HNEW = H*R
        ENDIF
        H=HNEW
      ELSE

%       If IPHASE = 0, increase order by one and multiply stepsize by
%       factor two.
        K = KP1
        HNEW = H*2.0D0
        H = HNEW
      ENDIF

%     Update differences for next step.
      IF (KOLD.NE.IWM(OMXORD)) CALL DCOPY(NEQ,ERRV,1,PHI(1,KP2),1)

%     Set PHI(*,KP1)=PHI(*,KP1)+ERRV(*).
      CALL DAXPY(NEQ,ONE,ERRV,1,PHI(1,KP1),1)

%     FOR J=KP1-1,1,-1 DO
%       SET PHI(*,J)=PHI(*,J)+PHI(*,J+1).
      DO 560 J=KP1-1,1,-1
         CALL DAXPY(NEQ,ONE,PHI(1,J+1),1,PHI(1,J),1)
 560  CONTINUE
      RETURN

C----------------------------------------------------------------------
C     Block 6
C     The step is unsuccessful. Restore X, PSI, PHI
C     determine appropriate stepsize for
C     continuing the integration, or exit with
C     an error flag if there have been many
C     failures.
C----------------------------------------------------------------------
 600  IPHASE = 1
C
C     restore T, PHI, PSI.
      T=TOLD
      if(KP1 >= NSP1)
C
C       FOR J=NSP1,KP1 DO
C         SET PHI(*,J)=PHI(*,J)/BETA(J).
        for J=NSP1:KP1
           TEMP1 = ONE/BETA(J);
           CALL DSCAL(NEQ,TEMP1,PHI(1,J),1)
        end
      end
      DO 620 I=2,KP1
        PSI(I-1)=PSI(I)-H
 620  CONTINUE

C     Test whether failure is due to problems with the
C     linear or least sqaure system, or the error test.
      IF((IER .EQ. 0).AND.(IRED.EQ.0))GO TO 660
      IWM(OCTF)=IWM(OCTF)+1
      IF(IER.EQ.0) GOTO 650

C     The lineare system is (maybe) singular. We got
C     an error in the factorization or solution subroutine
C     Reduce the stepsize by a factor of 4. If
C     this happens three times in a row on
C     the same step, return with an error flag.
      NSF=NSF+1
      R = 0.25D0
      H=H*R
      IF (NSF .LT. 3 .AND. ABS(H) .GE. HMIN) GO TO 690
      IDID=-8
      CALL DBDTRP(T,T,NEQ,K,PHI,LDPHI,PSI,X,XPRIME)
      RETURN
C
C     The computation of the solution of the system failed
C     for a reason other than a singular system. If IRED=-2,
C     then return. Otherwise, if IRED=-1 reduce the stepsize
C     and try again, unless too many failures have occured.
C     In all other cases we must return control to the calling
C     program.
 650  CONTINUE
      IF (IRED .EQ. -1) THEN
         NCF = NCF + 1
         R = 0.25D0
         H = H*R
         IF (NCF .LT. 10 .AND. ABS(H) .GE. HMIN) GOTO 690
         IDID = -10
      ELSEIF (IRED .EQ. -2) THEN
         IDID = -21
      ELSE
         IDID = IRED
          IF (IRED .EQ. -26) THEN
            M = MQ
            ID = IDQ
            IA = IAQ
            IU = IUQ
         ENDIF
      ENDIF
      CALL DBDTRP(T,T,NEQ,K,PHI,LDPHI,PSI,X,XPRIME)
      RETURN;
C
C     The cause of the failure was the error estimate
C     exceeding the tolerance.
 660  NEF=NEF+1;
      IWM(OETF)=IWM(OETF)+1
C
C     On first error test failure, keep current order or lower
C     order by one.  Compute new stepsize based on differences
C     of the solution.
      IF (NEF .EQ. 1) THEN
        K = KNEW
        TEMP2 = K + 1
        R = 0.90D0*(2.0D0*EST+0.0001D0)**(-ONE/TEMP2)
        R = MAX(0.25D0,MIN(0.9D0,R))
        H = H*R
C
C     On second error test failure, use the current order or
C     decrease order by one.  Reduce the stepsize by a factor of four.
      ELSEIF (NEF .EQ. 2) THEN
        K = KNEW
        H = 0.25D0*H
        
C     On third and subsequent error test failures, set the order to
C     one and reduce the stepsize by a factor of four.
      ELSE
        K = 1
        H = 0.25D0*H
      ENDIF

C     Check if stepsize is to small.
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
C
C     For all crashes, restore X to its last value,
C     interpolate to find XPRIME at last T, and return.
 675  CONTINUE
      CALL DBDTRP(T,T,NEQ,K,PHI,LDPHI,PSI,X,XPRIME)
      RETURN
C
C     go back and try this step again.
 690  GO TO 200
%}