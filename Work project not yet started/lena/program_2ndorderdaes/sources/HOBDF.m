% ONE STEP !!!!!!!!!
% N number of equations
% T initial point (overwritten)
% TOUT
% X initial values (overwritten)
% RTOL, ATOL error tolerances scaler or vectors
% METHOD integration method 1 =BDF, 2 =RK, 3=GLM
%
function [T,X,XPRIME,H,RWORK,IWORK,IERR]=HOBDF(F,dF,ddF,NEQ,T,TOUT,X,XPRIME,RTOL,ATOL,METHOD,INFO,LRW,LIW,RWORK)

%Parameters
% OM=1; 
% OID=2; 
% OIA=3; 
% OIU=4;
% OIWM=1; 

% Indices von IWORK
OMXORD=3; 
OMXINDX=4; 
% OPHASE=6; 
OK=7; 
% OKOLD=8; 
% ONS=9; 
ONSTL=10; 
ONST=11; 
% ONEV=12; 
% ONFA=13; 
OETF=14; 
OCTF=15; 
% ONSING=16;
% OREQU=17; 
% OCEQU=18; 
% OIREQ=19; 
 ONMAX=20;

%Indices von RWORK
OTSTOP=1; 
OHMAX=2; 
OH=3; 
OTN=4; 
% OCJ=5; 
% OCJOLD=6;
% OHOLD=7; 
% OERRACC=8; 
OROUND=9; 
OSAFE=11;
OFACL=12; 
OFACR=13; 
OQUOT1=14; 
OQUOT2=15; 
% OE=21;


if( INFO(1)==0)
    %----------------------------------------------------------------------
    %     This block is executed for the initial call only.
    %     It contains checking of inputs and initializations.
    %----------------------------------------------------------------------

    %     First check INFO array to make sure all elements of INFO
    %     are either zero or one.
    for i=2:20
        if (INFO(i) ~= 0 && INFO(i)~=1 )
            IERR = -101;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
    end

    %     Check if NEQ is greater than zero.
    if (NEQ <= 0 )
        IERR = -102;
        if (INFO(1) == -1)
            IERR = -999;
        end
        INFO(1) = -1;
        return;
    end

    %     Check and compute maximum order.
    MXORD = 5;
    if (INFO(9)~=0)
        MXORD = IWORK(OMXORD);
        if (MXORD<1 || MXORD>5)
            IERR = -103;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
    end
    IWORK(OMXORD) = MXORD;

    %   Check and compute maximum index.
    MXINDX=3;
    if (INFO(10)~=0)
        MXINDX=IWORK(OMXINDX);
        if (MXINDX<0)
            IERR = -104;
            if(INFO(1)==-1)
                IERR = -999;
            end

            INFO(1)=-1;
            return;
        end
    end
    IWORK(OMXINDX)=MXINDX;
    
    %     Compute LENRW and LENIW.
    if (METHOD == 1)
        LENRW=50+(21+IWORK(OMXORD)+14*IWORK(OMXINDX))*NEQ + (6+IWORK(OMXINDX)+6*(IWORK(OMXINDX)+1)^2)*NEQ*NEQ;
        LENIW=20+2*NEQ;
    else
        LENRW=20 + (15+14*(IWORK(OMXINDX)+1))*NEQ+(14+(IWORK(OMXINDX)+1)+6*(IWORK(OMXINDX)+1)^2)*NEQ*NEQ;
        LENIW=20+6*NEQ;
    end
    
    %     Check lengths of RWORK and IWORK.
    if (LRW<LENRW)
        IWORK(LRW) = LENRW;
        IERR = -105;
        if(INFO(1)==-1)
            IERR = -999;
        end

        INFO(1)=-1;
        return;

    end
    if (LIW<LENIW)
        IWORK(LIW) = LENIW;
        IERR = -106;
        if(INFO(1)==-1)
            IERR = -999;
        end

        INFO(1)=-1;
        return;
    end


    %     Check HMAX.
    if (INFO(7)~=0)
        HMAX = RWORK(OHMAX);
        if (HMAX<=0.0)
            IERR = -111;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
    end

    %     Check NMAX , the maximal number of steps.
    if (INFO(6)==0)
        NMAX = 10000;
    else
        NMAX = IWORK(ONMAX);
        if (NMAX<=0)
            IERR = -120;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
    end
    IWORK(ONMAX) = NMAX;

    %     Check the safety factor and parameters in step size prediction.
    if (INFO(16)==0)      % INFO(16) ????????
        SAFE=0.90;
        FACL=5.00;
        FACR=0.1250;
        RWORK(OSAFE)=SAFE;
        RWORK(OFACL)=FACL;
        RWORK(OFACR)=FACR;
    else
        SAFE=RWORK(OSAFE);
        if (SAFE<=0.001 || SAFE>= 1.00)
            IERR = -121;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
        FACL=RWORK(OFACL);
        FACR=RWORK(OFACR);
        if (FACL < 1.0 || FACR > 1.0)
            IERR = -122;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
    end

    %     QUOT1 and QUOT2: If QUOT1 < HNEW/HOLD < QUOT2, stepsize = const.
    if (INFO(17)==0)    % INFO(17) ???
        QUOT1=1.00;
        QUOT2=1.20;
        RWORK(OQUOT1)=QUOT1;
        RWORK(OQUOT2)=QUOT2;
    else
        QUOT1=RWORK(OQUOT1);
        QUOT2=RWORK(OQUOT2);
        if (QUOT1>1.00 || QUOT2<1.00)
            IERR = -123;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
    end


    %     Initialize counters
    IWORK(ONST)  = 0;
    %IWORK(ONRE)  = 0;
    %IWORK(ONJE)  = 0;
    IWORK(OETF)  = 0;
    IWORK(OCTF)  = 0;
    IWORK(ONSTL) = 0;

     %IWORK(ONEV)=0;
     %IWORK(ONFA)=0;

    IWORK(OK)    = 1;
    IERR = 1;
else
    % ----------------------------------------------------------------------
    %      This block is for continuation calls
    %      only. Here we check INFO(1), and if the
    %      last step was interrupted we check whether
    %      appropriate action was taken.
    % ----------------------------------------------------------------------

    if (INFO(1)>=1)
        IWORK(ONSTL) = IWORK(ONST);
    else
        if (INFO(1)~=-1)
            IERR = -101;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end

        %      If we are here, the last step was interrupted
        %      by an error condition from DBDFST, and
        %      appropriate action was not taken. This
        %      is a fatal error.
        print(' THE LAST STEP TERMINATED WITH A NEGATIVE VALUE OF', ' IERR =',IERR);
        print(' AND NO APPROPRIATE ACTION WAS TAKEN.');
        print(' RUN TERMINATED');
        IERR = -998;
        stop; %return;
    end
end


% ----------------------------------------------------------------------
%      This block is executed on all calls.
%      The error tolerance parameters are
%      checked, and the work array pointers
%      are set.
% ----------------------------------------------------------------------

%      Check if E(t) and A(t) are constant.
if (INFO(5)==0)
    MCONST = false;
else
    MCONST = true;
end

%     Check whether predictive step size control shall be used.
if (INFO(15)==0)  % INFO(15)=?
    PRED = true;
else
    PRED = false
end

%      Check  RTOL and ATOL.
NZFLG = 0;
if (INFO(2)==0)
    RTOLI = RTOL(1);
    ATOLI = ATOL(1);
    if (RTOLI> 0.0 || ATOLI>0.0)
        NZFLG=1;
    end
    if (RTOLI<0.0)
        IERR = -107;
        if (INFO(1) == -1)
            IERR = -999;
        end
        INFO(1) = -1;
        return;
    end
    if (ATOLI<0.0)
        IERR = -108;
        if (INFO(1) == -1)
            IERR = -999;
        end
        INFO(1) = -1;
        return;
    end
else
    for i=1:NEQ
        RTOLI = RTOL(i);
        ATOLI = ATOL(i);
        if (RTOLI>0.0 || ATOLI>0.0)
            NZFLG=1;
        end
        if (RTOLI<0.0)
            IERR = -107;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
        if (ATOLI<0.0)
            IERR = -108;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
    end
end

if (NZFLG==0)
    IERR = -109;
    if (INFO(1) == -1)
        IERR = -999;
    end
    INFO(1) = -1;
    return;
end

%      Set up leading dimensions for RWORK matrices.
% if (METHOD == 1)
%     LDPHI   = NEQ;
%     LDW     = NEQ;
% else
%     N3      = 3*NEQ;
%     LDCONT  = NEQ;
%     LDW     = N3;
%     LDE0    = NEQ;
%     LDA0    = NEQ;
% end
% LDE     = NEQ;
% LDA     = NEQ;
% LDZ1    = NEQ;
% LDZ2Q   = (IWORK(OMXINDX)+1)*NEQ;
% LDEQ    = (IWORK(OMXINDX)+1)*NEQ;
% LDAQ    = (IWORK(OMXINDX)+1)*NEQ;
% LDAH    = (IWORK(OMXINDX)+1)*NEQ;
% 
% LJACQ   = N + ID;
% LZ2Q    = N;
% LZ1     = NEQ;
% LDPHI   = NEQ;

%     Set up RWORK storage. IWORK is fixed in DATA SEGMENT.
% OA      = OE + LDE*LDE;
% OF      = OA + LDA*LDA;
% OEQ     = OF + NEQ;
% OAQ     = OEQ + LDEQ*LDEQ;
% OFQ     = OAQ + LDAQ*NEQ;
% OZ1     = OFQ + (IWORK(OMXINDX)+1)*NEQ;
% OZ2Q    = OZ1 + LDZ1*LDZ1;
% OAH     = OZ2Q + LDZ2Q*LDZ2Q;
% OD      = OAH + LDAH*LDAH;
% OV      = OD + (IWORK(OMXINDX)+1)*NEQ;
% OWT     = OV +  (IWORK(OMXINDX)+1)*NEQ;
% if (METHOD == 1)
%     OALPHA  = OWT + NEQ;
%     OBETA   = OALPHA + 6;
%     OGAMMA  = OBETA + 6;
%     OPSI    = OGAMMA + 6;
%     OSIGMA  = OPSI + 6;
%     ODELTA  = OSIGMA + 6;
%     OERRV   = ODELTA + NEQ;
%     OPHI    = OERRV + NEQ;
%     ORS     = OPHI + (IWORK(OMXORD)+1)*LDPHI;
%     OCS     = ORS + NEQ;
%     OW      = OCS + NEQ;
% else
%     OZR1    = OWT + NEQ;
%     OZR2    = OZR1 + NEQ;
%     OZR3    = OZR2 + NEQ;
%     OYY     = OZR3 + NEQ;
%     OB      = OYY + N3;
%     OCONT   = OB + N3;
%     OE0     = OCONT + LDCONT*4;
%     OA0     = OE0 + LDE0^2;
%     OW      = OA0 + LDA0^2;
% end
% OWORK   = OW + LDW*LDW;
% LWORK   = LRW - OWORK;


if(INFO(1)==0)
    % ----------------------------------------------------------------------
    %      This block is executed on the initial call  only
    %      Set the initial step size, and the error weight vector XW, and PHI.
    %      Compute initial XPRIME and X, if necessary.
    % ---------------------------------------------------------------------

    TN = T;
    IERR = 1;

    %     Compute unit roundoff.
    UROUND = eps; %DLAMCH('P');
    RWORK(OROUND) = UROUND;
    if (INFO(1)==0)
        if (TOUT == T)
            INFO(1) = 2;
            IERR = 4;

            RWORK(OTN)=TN;
            RWORK(OH)=H;
            RETURN;
        end
    end

    %     Check to see that TOUT is different from T.
    if (TOUT == T)
        IERR = -119;
        if (INFO(1) == -1)
            IERR = -999;
        end
        INFO(1) = -1;
        return;
    end

    %     Compute HMIN.
    HMIN = 4.0*UROUND*max(abs(T),abs(TOUT));

    %     Check initial interval to see that it is long enough.
    TDIST = abs(TOUT - T);
    if (TDIST < HMIN)
        IERR = -115;
        if (INFO(1) == -1)
            IERR = -999;
        end
        INFO(1) = -1;
        return;
    end

    %     Set error weight vector WT.
    RTOLI = RTOL(1);
    ATOLI = ATOL(1);
    for i=1:NEQ
        if (INFO(2) ~=0)
            RTOLI = RTOL(i);
            ATOLI = ATOL(i);
        end
       % RWORK(OWT+i-1) = RTOLI*abs(X(i)) + ATOLI;
       WT(i)=RTOLI*abs(X(i)) + ATOLI;
       % if (RWORK(OWT+i-1)<=0.0)
       if (WT(i)<=0.0)
            IERR = -114;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
    end

    %     Invert the WT vector to minimize number of division operations.
    for i = 1:NEQ
        %RWORK(OWT+i-1) = 1.0/RWORK(OWT+i-1);
        WT(i) = 1.0/WT(i);
    end

    %     Check H0, if this was input.
    if (INFO(8) ~= 0)
        H0 = RWORK(OH);
        if ((TOUT - T)*H0 < 0.0)
            IERR = -112;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
        if (H0 == 0.0)
            IERR = -113;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
    else
        %     Compute initial stepsize, to be used by DBDFST.
        H0 = 0.001D0*TDIST;
        %XPNORM = DGENRM(NEQ,X(NEQ+1),RWORK(OWT));
        XPNORM = DGENRM(NEQ,X(NEQ+1),WT);
        if (XPNORM > 0.5/H0)
            H0 = 0.5/XPNORM;
        end
        H0 = abs(H0)*sign(TOUT-T);
    end

    %     Adjust H0 if necessary to meet HMAX bound.
    if (INFO(7) ~= 0)
        RH = abs(H0)/RWORK(OHMAX);
        if (RH > 1.0)
            H0 = H0/RH;
        end
    end

    %     Compute TSTOP, if applicable.
    if (INFO(4) ~= 0)
        TSTOP = RWORK(OTSTOP);
        if ((TSTOP - T)*H0 < 0.0)
            IERR = -116;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
        if ((T + H0 - TSTOP)*H0 > 0.0)
            H0 = TSTOP - T;
        end
        if ((TSTOP - TOUT)*H0 < 0.0)
            IERR = -110;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
    end

    %     Load H with H0. store H in RWORK(OH).
    H = H0;
    RWORK(OH) = H;

    
    if (METHOD == 1)
        %DCOPY(NEQ,X,1,RWORK(OPHI),1)
        PHI = X;
        %ITEMP = OPHI + NEQ;
        for i = 1:NEQ
            %RWORK(ITEMP + i - 1) = H*XPRIME(i);
            %RWORK(ITEMP + N + i - 1) = 0.0;
            PHI(i) = H*XPRIME(i);           %??????
            PHIPRIME(i) = 0.0;              % ?????
        end
    end

    %360
    if (IERR<0)
        INFO(1) = -1;
        T = TN;
        RWORK(OTN) = TN;
        RWORK(OH)  = H;
        return;
    end
end



if(INFO(1) == 1)
    % -------------------------------------------------------
    %      This block is for continuation calls only. Its
    %      purpose is to check stop conditions before taking a step.
    %      Adjust H if necessary to meet HMAX bound.
    % -------------------------------------------------------
    UROUND=RWORK(OROUND);
    DONE = false;
    TN = RWORK(OTN);
    H = RWORK(OH);
    if (INFO(7) ~= 0)
        RH = abs(H)/RWORK(OHMAX);
        if (RH > 1.0)
            H = H/RH;
        end
    end

    if (T == TOUT)
        IERR = -119;
        if (INFO(1) == -1)
            IERR = -999;
        end
        INFO(1) = -1;
        return;
    end
    if ((T - TOUT)*H > 0.0)
        IERR = -112;
        if (INFO(1) == -1)
            IERR = -999;
        end
        INFO(1) = -1;
        return;
    end
    if (INFO(4) == 1)
        TSTOP = RWORK(OTSTOP);
        if ((TN-TSTOP)*H > 0.0)
            IERR = -116;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
        if ((TSTOP-TOUT)*H < 0.0)
            IERR = -110;
            if (INFO(1) == -1)
                IERR = -999;
            end
            INFO(1) = -1;
            return;
        end
    end

    if (METHOD == 1)
        %        Backward-differentiation branch.
        if (INFO(4) ~= 1)
            if (INFO(3) ~= 1)
                if ((TN-TOUT)*H < 0.0)
                    if (DONE)
                        % All successful returns from DGENDA are made from this block.
                        RWORK(OTN) = TN;
                        RWORK(OH)  = H;
                        return;
                    end
                end
               % DBDTRP (TN, TOUT,NEQ,IWORK(OKOLD), RWORK(OPHI),LDPHI,RWORK(OPSI),X,XPRIME);
                DBDTRP (TN, TOUT,NEQ,IWORK(OKOLD), PHI,LDPHI,RWORK(OPSI),X,XPRIME);
                T=TOUT;
                IERR = 3;
            else
                if ((TN-T)*H <= 0.0)
                    if (DONE)
                        % All successful returns from DGENDA are made from this block.
                        RWORK(OTN) = TN;
                        RWORK(OH)  = H;
                        return;
                    end
                end
                if ((TN - TOUT)*H <= 0.0)
                    %DBDTRP (TN, TN, N, IWORK(OKOLD), RWORK(OPHI), LDPHI,RWORK(OPSI), X, XPRIME);
                    DBDTRP (TN, TN, N, IWORK(OKOLD), PHI, LDPHI,RWORK(OPSI), X, XPRIME);
                    T = TN;
                    IERR = 1;
                else
                    %DBDTRP (TN, TOUT, N, IWORK(OKOLD), RWORK(OPHI), LDPHI, RWORK(OPSI), X, XPRIME);
                    DBDTRP (TN, TOUT, N, IWORK(OKOLD), PHI, LDPHI, RWORK(OPSI), X, XPRIME);
                    T = TOUT;
                    IERR = 3;
                end
            end

        elseif (INFO(3) ~= 1)
            if ((TN-TOUT)*H<0.0)
                print('hier');
                error('weiter?');
                %GOTO 410
            end
            %DBDTRP (TN, TOUT, N, IWORK(OKOLD), RWORK(OPHI), LDPHI,    RWORK(OPSI), X, XPRIME);
            DBDTRP (TN, TOUT, N, IWORK(OKOLD), PHI, LDPHI,    RWORK(OPSI), X, XPRIME);
            T=TOUT;
            IERR = 3;
        else
            if ((TN-T)*H <= 0.0)
                print('hier');
                error('weiter?');
                %GOTO 410
            end
            if ((TN - TOUT)*H <= 0.0)
                %DBDTRP(TN, TN, N, IWORK(OKOLD), RWORK(OPHI), LDPHI, RWORK(OPSI), X, XPRIME)
                DBDTRP(TN, TN, N, IWORK(OKOLD), PHI, LDPHI, RWORK(OPSI), X, XPRIME)
                T = TN;
                IERR = 1;
            else
                %DBDTRP(TN, TOUT, N, IWORK(OKOLD), RWORK(OPHI), LDPHI,RWORK(OPSI), X, XPRIME)
                DBDTRP(TN, TOUT, N, IWORK(OKOLD), PHI, LDPHI,RWORK(OPSI), X, XPRIME)
                T = TOUT;
                IERR = 3;
            end
        end
        DONE = true;
        if (DONE)
            %      All successful returns from DGENDA are made from this block.
            RWORK(OTN) = TN;
            RWORK(OH)  = H;
            return;
        end

        %      Check whether we are within roundoff of tstop.
        if (abs(TN-TSTOP) <= 100.0*UROUND*(abs(TN)+abs(H)))
            %DBDTRP (TN, TSTOP, N, IWORK(OKOLD), RWORK(OPHI), LDPHI, RWORK(OPSI), X, XPRIME);
            DBDTRP (TN, TSTOP, N, IWORK(OKOLD), PHI, LDPHI, RWORK(OPSI), X, XPRIME);
            IERR = 2;
            T = TSTOP;
            DONE = true;
            if (DONE)
                %      All successful returns from DGENDA are made from this block.
                RWORK(OTN) = TN;
                RWORK(OH)  = H;
                return;
            end
        else
            TNEXT = TN + H;
            if ((TNEXT-TSTOP)*H <= 0.0)
                if (DONE)
                    %      All successful returns from DGENDA are made from this block.
                    RWORK(OTN) = TN;
                    RWORK(OH)  = H;
                    return;
                end
            end
            H = TSTOP - TN;
            RWORK(OH) = H;
        end
    end

    % end if BDF  ??
    %420
    if (DONE)
        RWORK(OTN)=TN;
        RWORK(OH)=H;
        RETURN;
    end
end


% -------------------------------------------------------
%      The next block contains the call to the one-step integrator DBDFST.
%      This is a looping point for the integration steps. Check for too many steps.
%      Update WT. Check for too much accuracy requested. Compute minimum stepsize.
% -------------------------------------------------------

%  500 
DONE = false;

%      Check for too many steps.
if ((IWORK(ONST)-IWORK(ONSTL)) >= IWORK(ONMAX))
    IERR = -1;
if (IERR < 0)
    INFO(1) = -1;
    T = TN;
    RWORK(OTN) = TN;
    RWORK(OH)  = H;
    return;
end
        
end

%      Update WT.
RTOLI = RTOL(1);
ATOLI = ATOL(1);
for i=1:NEQ
    if (INFO(2) ~=0)
        RTOLI = RTOL(i);
        ATOLI = ATOL(i);
    end
    %RWORK(OWT+i-1) = RTOLI*abs(RWORK(OMXORD+i-1)) + ATOLI;
    WT(i) = RTOLI*abs(RWORK(OMXORD+i-1)) + ATOLI;
    %if (RWORK(i+OWT-1) <= 0.0)
    if (WT(i) <= 0.0)
        IERR = -3;
if (IERR < 0)
    INFO(1) = -1;
    T = TN;
    RWORK(OTN) = TN;
    RWORK(OH)  = H;
    return;
end
    end
end

%      Invert the WT vector to minimize number of division operations.
for i=1:NEQ
    %RWORK(OWT+i-1) = 1.0/RWORK(OWT+i-1);
    WT(i) = 1.0/WT(i);
end

%      Test for too much accuracy requested.
%R = DGENRM (NEQ, RWORK(OMXORD), RWORK(OWT)) * 100.0 * UROUND;
R = DGENRM (NEQ, RWORK(OMXORD), WT) * 100.0 * UROUND;
if (R>1.0)
    %     %        Multiply RTOL and ATOL by R and return.
    if (INFO(2)~=1)
        RTOL(1) = R*RTOL(1);
        ATOL(1) = R*ATOL(1);
    else
        DSCAL(NEQ, R, RTOL, 1);
        DSCAL(NEQ, R, ATOL, 1);
    end
    IERR = -2;
    if (IERR < 0)
        INFO(1) = -1;
        T = TN;
        RWORK(OTN) = TN;
        RWORK(OH)  = H;
        return;
    end
end

%      Compute minimum stepsize.
HMIN = 4.0 * UROUND * max(abs(TN), abs(TOUT));

%      Test H vs. HMAX.
if (INFO(7) ~= 0)
    RH = abs(H)/RWORK(OHMAX);
    if (RH > 1.0)
        H = H/RH;
    end
end

%     Runge-kutta integrator should not step past TOUT.
if (METHOD == 2 && (TN+H*1.0001-TOUT)*H >= 0.0)
    HO=H;
    H=TOUT-TN;
    DONE=true;
end

if (METHOD == 1)
%    DBDSTP(EDIF, ADIF, FDIF, NEQ, TN, H, RWORK(OHOLD), HMIN,
%    RWORK(OCJ),RWORK(OCJOLD), IWORK(OK), IWORK(OKOLD), CVAL(OM), CVAL(OID), CVAL(OIA), CVAL(OIU),
%   IWORK(OIREQ), IWORK(OMXINDX), IWORK(ONS),IWORK(OPHASE), IWORK(OREQU), IWORK(OCEQU), INFO(1),
%    X, XPRIME, RWORK(OE), LDE, RWORK(OA), LDA, RWORK(OZ1), LDZ1, RWORK(OZ2Q),
%    LDZ2Q, RWORK(OWT), RWORK(OW), LDW, RWORK(OPHI), LDPHI, RWORK(OALPHA),
%   RWORK(OBETA), RWORK(OGAMMA), RWORK(OPSI), RWORK(OSIGMA), RWORK(ODELTA), RWORK(OERRV),
%    RWORK(ORS), RWORK(OCS), RWORK(OF), RWORK(OEQ), LDEQ, RWORK(OAQ), LDAQ, RWORK(OFQ), RWORK(OAH),
%    LDAH, IPAR, RPAR, IWORK(OIWM), RWORK(OWORK), LWORK, MCONST, IWARN, IERR);
else
%    DRKSTP(EDIF, ADIF, FDIF, NEQ, N3, TN, H, RWORK(OHOLD),HMIN, RWORK(OERRACC), RWORK(OROUND), RWORK(OSAFE),
%    RWORK(OFACL), RWORK(OFACR), RWORK(OQUOT1),RWORK(OQUOT2), CVAL(OM), CVAL(OID), CVAL(OIA),
%   CVAL(OIU), IWORK(OIREQ), IWORK(OMXINDX),IWORK(ONSING), INFO(1), DONE, X, XPRIME,
%    RWORK(OE0), LDE0, RWORK(OA0), LDA0, RWORK(OE), LDE,RWORK(OA), LDA, RWORK(OEQ), LDEQ, RWORK(OAQ), LDAQ,
%    RWORK(OZ1), LDZ1, RWORK(OZ2Q), LDZ2Q, RWORK(OAH),LDAH, RWORK(OW), LDW, RWORK(OCONT), LDCONT,
%   RWORK(OWT), RWORK(OZR1), RWORK(OZR2), RWORK(OZR3),RWORK(OYY), RWORK(OB), RWORK(OF), RWORK(OFQ), IPAR,
%    RPAR, IWORK(OIWM), RWORK(OWORK), LWORK, PRED,MCONST, IWARN, IERR);
end


if (IERR < 0)
    INFO(1) = -1;
    T = TN;
    RWORK(OTN) = TN;
    RWORK(OH)  = H;
    return;
end

% %     Consistency check for inhomogeneity
% MY_DINSYS(NQ, MQ, 0, RWORK(OXQ), RWORK(OFQ), IERR, TN, M, IM,  FDIF, IPAR, RPAR);
%
%
% DCOPY (NQ, X, 1, RWORK(OXQOLD), 1);
% DCOPY (NQ, RWORK(OXQ), 1, X, 1);
% if (INFO(16) ~= 0)
%     if (NEWTON == 0)
%         IERR = -126;
%         if (INFO(1) == -1)
%             IERR = -999;
%         end
%         INFO(1) = -1;
%         return;
%     end
%     if ((INFO(17)==0 && (IWORK(ONST)-IWORK(ONSTL))==1)|| INFO(17)==1)
%         IDID = 0;
%         MY_DCSIND (SCALE, DFDIF, USCAL, M, N, TN, IWORK(OMXINDX),  X, RWORK(OCOND), IPAR, RPAR, SCALC, SCALR, IM,   ID, IA, IU, RWORK(OJACQ), LJACQ,  RWORK(ORWORK), LRWORK, IERR)
%         if (IERR<0)
%             INFO(1) = -1;
%             T = TN;
%             RWORK(OTN) = TN;
%             RWORK(OH)  = H;
%             return;
%         end
%     end

% --------------------------------------------------------
%      This block handles the case of a successful return
%      from DBDFST (IERR=1).
%      Test for stop conditions.
% --------------------------------------------------------

%     Runge-Kutta integrator reached exactly TOUT.
if (DONE)
    H=HO;
    T=TN;
    IERR=3;
    RWORK(OTN)=TN;
    RWORK(OH)=H;
    RETURN;
end

%     Backward-differentiation branch.
if (METHOD == 1)
    if (INFO(4) == 0)
        if (INFO(3) == 0)
            if ((TN-TOUT)*H < 0.0)
                %GOTO 500
            end
%            DBDTRP (TN, TOUT, NEQ, IWORK(OKOLD), RWORK(OPHI), LDPHI,  RWORK(OPSI), X, XPRIME);
            IERR = 3;
            T = TOUT;
        else
            if ((TN-TOUT)*H < 0.0)
                T = TN;
                IERR = 1;
            else
  %              DBDTRP (TN, TOUT, NEQ, IWORK(OKOLD), RWORK(OPHI), LDPHI,   RWORK(OPSI), X, XPRIME);
                IERR = 3;
                T = TOUT;
            end
        end
    else
        if (INFO(3) == 0)
            if ((TN-TOUT)*H >= 0.0)
    %            DBDTRP (TN, TOUT, NEQ, IWORK(OKOLD), RWORK(OPHI), LDPHI,  RWORK(OPSI), X, XPRIME);
                T = TOUT;
                IERR = 3;
            else
                if (abs(TN-TSTOP) > 100.0*UROUND*(abs(TN)+abs(H)))
                    TNEXT = TN+H;
                    if ((TNEXT-TSTOP)*H <= 0.0)
                        %GO TO 500
                    end
                    H = TSTOP-TN;
                    %GO TO 500
                else
        %            DBDTRP (TN, TSTOP, NEQ, IWORK(OKOLD), RWORK(OPHI), LDPHI,  RWORK(OPSI), X, XPRIME);
                    IERR = 2;
                    T = TSTOP;
                end
            end

        else
            if ((TN-TOUT)*H < 0.0)
                if (abs(TN-TSTOP) > 100.0*UROUND*(abs(TN)+abs(H)))
                    T = TN;
                    IERR = 1;
                else
                    DBDTRP (TN, TSTOP, NEQ, IWORK(OKOLD), RWORK(OPHI), LDPHI,  RWORK(OPSI), X, XPRIME);
                    IERR = 2;
                    T = TSTOP;
                end
            else
                DBDTRP (TN, TOUT, NEQ, IWORK(OKOLD), RWORK(OPHI), LDPHI, RWORK(OPSI), X, XPRIME);
                T = TOUT;
                IERR = 3;
            end
        end
    end
else  %     Runge-Kutta branch
    if (INFO(3) == 0)
        GOTO 500
    else
        T=TN;
        IERR=1;
    end
end

%--------------------------------------------------------
%     All successful returns from DGELDA are made from
%     this block.
%--------------------------------------------------------

%580   CONTINUE
RWORK(OTN)=TN;
RWORK(OH)=H;
return;