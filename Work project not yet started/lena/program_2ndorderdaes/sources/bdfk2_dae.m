%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BDFk-Verfahren, variable-step variable-order BDF
%
% Loesen der differentiell-algebraischen Gleichung 2. Ordung der Form
% M y''(t)=D y'(t) + S y(t) + f(t,y,y')
%
% mit BDF der Ordnung Ordnung k
%
% CALL  : [t,h,y,y_est,y_interpol]=bdfk2dae(funct,y0,y1,t0,a,N,K,RTOL,ATOL,H0)
%
% INPUT : funct - Name der aufzurufenden Routine, die die Matrizen M, D, S und Vektor f bestimmt
%                 (Text-String)
%         y0,y1 - initial value
%         t0    - initial time         t\in[t0,t0+a]
%         a     - length of interval   t\in[t0,t0+a]
%         N     - number of steps
%         K     - maximal order
%
% OUTPUT: t    - time steps t=[t0,t0+h,...,t0+(N-1)*h,t0+N*h=t0+a]
%         h    - stepsizes
%         y    - approximation of solution at t(i)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,h,y]=bdfk2_dae(funct,funct_d1,funct_d2,funct_d3,nL,y0,y1,t0,a,N,K,RTOL,ATOL,H0,step)

% check input arguments
maxOrder=3;
if(K>maxOrder || K<1)
     error('wrong order');
end

if(a==0)
    error('wrong interval');
end

TIME=t0:a/N:t0+a;
t(1)=t0;

TOUT=t0+a;

% initial stepsize
if(H0==0)
    H0 = sign(TOUT-t0)*min(10^(-3)*abs(TOUT-t0),0.5*myNorm(y1,RTOL,ATOL)^(-1));
end

h(1)=H0(1);

[n,m]=size(y0);
[n1,m1]=size(y1);

y(1:n,1)=y0;
y(n+1:n+n1,1)=y1;
y_interpol(:,1)=y(:,1);

iNoRej=0;nef=0;ncf=0;nsf=0;iPhase=0;newOrder = 0;

% integration
i=1;I=1;

while(i<=2000 && t(i)<=TOUT)
    i
    t(i+1) = t(i)+h(i);
    h(i+1) = h(i);

    %order selection -> k
    if( i <= 11)
        if(K==1)
            if(i==1)
                KK(1)=1;
                KK(2)=2;
                KK(3)=0;
            end
            if(i>1)
                KK(1)=1;
                KK(2)=1;
                KK(3)=0;
            end
        end
        if(K==2)
            if(i==1)
                KK(1)=1;
                KK(2)=2;
                KK(3)=0;
            end
            if(i==2)
                KK(1)=2;
                KK(2)=1;
                KK(3)=3;
            end
            if(i==3)
                KK(1)=2;
                KK(2)=2;
                KK(3)=1;
            end
            if(i>3)
                KK(1)=2;
                KK(2)=2;
                KK(3)=2;
            end
        end
    else %( i> 4)
        if(newOrder==1) % increase order
            if( KK(1) == maxOrder) % k=k+1
                temp = KK(1);
            else
                temp = KK(1)+1;
            end
            for zz = temp+1:-1:2
                KK(zz) = KK(zz-1);
            end
            KK(1) = temp;

        elseif(newOrder==2) % decrease order
            if( KK(1) == 1)
                KK(1) = KK(1);
            else
                KK(1) = KK(1)-1;
            end
            for zz = KK(1)+1:-1:2
                KK(zz) = KK(zz-1) ;
            end
        elseif(newOrder==3) % order constant
            KK(1) = KK(1);
            for zz=KK(1)+1:-1:2
                KK(zz) = KK(zz-1);
            end
        else
            error('order selection');
        end
    end

    
    if( i<4 )
        if(i==1) 
            T=[t(2),t(1),t(1)-h(1),t(1)-2*h(1),t(1)-3*h(1)];
        end
        if(i==2) 
            T=[t(3),t(2),t(1),t(1)-h(1),t(1)-2*h(1),t(1)-3*h(1)];
        end
        if(i==3) 
            T=[t(4),t(3),t(2),t(1),t(1)-h*(1)];
        end
    else
        time = KK(KK(1)+1)+KK(1)+1;
        T=t(i+1:-1:i- time +2);
    end


    [y(:,i+1),converged, IER ] = bdf2_step(funct,funct_d1,funct_d2,funct_d3,n,nL,KK,T,t,h,i,y,RTOL,ATOL);


    % Error Estimate - extrapolation
    if(i<=4)
        y_est(:,i+1)=extrapol(t,y,i-1);
        err2(:,i+1)=myNorm(y_est(1:n,i+1)-y(1:n,i+1),RTOL,ATOL);
    end
    if(i>4)
        y_est(:,i+1)=extrapol(t,y,3);
        err2(:,i+1)=myNorm(y_est(1:n,i+1)-y(1:n,i+1),RTOL,ATOL);
    end

    hmin = 4*eps*max(abs(t(i+1)),abs(TOUT));

    EST=err2(:,i+1);


    if( step ==1 ) % variable stepsize

        if( i>5)
            if( err2(:,i+1) <= 1.0)

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % step successful
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(iPhase==0)
                    newOrder = 1;
                    h(i+1)=2*h(i);
                    iPhase=0;
                    iNoRej=0;
                else
                    r=(2*EST+0.0001)^(-1/(KK(1)+1));
                    if( r < 2.0)
                        if( r > 1.0)
                            h(i+1)=h(i);
                        else
                            r = max(0.5,min(0.9,r));
                            h(i+1)=r*h(i);
                        end
                    else
                        h(i+1)=2*h(i);
                    end
                    iPhase=0;
                    newOrder = 3;
                end

                nef=0;
                ncf=0;
                nsf=0;

                % interpolation between meshpoints
                if(TIME(I+1)==t(i+1))
                    y_interpol(:,I+1)=y(:,i+1);
                    I=I+1;
                elseif(TIME(I+1)<t(i+1))
                    while(TIME(I+1)<t(i+1))
                        if(i>5)
                            nn=5;
                        else
                            nn=i;
                        end
                        y_interpol(:,I+1)=interpol(t,y,TIME(I+1),nn);
                        I=I+1;
                        if(I-1==N)
                            return
                        end
                    end
                end

            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % reject step, step unsucessful
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                iPhase=1;

                %restore y,t,y_est,err1
                y=y(:,1:i);
                t=t(1:i);
                y_est=y_est(:,1:i);
                err2=err2(:,1:i);

                i=i-1;
                iNoRej=iNoRej+1;

                %test whether failure was due to corrector iteration  or error test failure
                if(converged==0)   % Newton iteration converged
                    nef = nef+1;
                    if(nef > 1)
                        if( nef > 2)
                            newOrder=2;
                            r=0.25;
                            h(i+1)=r*h(i+1);
                            if(abs(h(i+1))<hmin)
                                error('minimal stepsize 1');
                            end
                        else % second error test failure
                            newOrder=3;
                            r=0.25;
                            h(i+1)=r*h(i+1);
                            if(abs(h(i+1))<hmin)
                                error('minimal stepsize 2');
                            end
                        end
                    else  % first error test failure
                        newOrder=3;
                        r=0.9*(2.0*EST+0.0001)^(-1/(KK(1)+1));
                        r=max(0.25,min(0.9,r));
                        h(i+1)=r*h(i+1);
                        if(abs(h(i+1))<hmin)
                            error('minimal stepsize 3');
                        end
                    end

                else % Newton iteration failed to converge
                    if(IER == 0)
                        ncf=ncf+1;
                        r=0.25;
                        h(i+1)=r*h(i+1);
                        if(abs(h(i+1))<hmin)
                            error('minimal stepsize 4');
                        end
                    else
                        nsf=nsf+1;
                        r=0.25;
                        h(i+1)=r*h(i+1);
                        if(abs(h(i+1))<hmin || nsf >3 )
                            error('minimal stepsize or nsf >3');
                        end

                    end

                    newOrder=3;
                end

            end
        else
            % interpolation between meshpoints
            if(TIME(I+1)==t(i+1))
                y_interpol(:,I+1)=y(:,i+1);
                I=I+1;
            elseif(TIME(I+1)<t(i+1))
                y_interpol(:,I+1)=interpol(t,y,TIME(I+1),i);
                I=I+1;
            end
        end

    else    % constant stepsize
        h(i+1) =h(i);


        if( i>10)
%            if( err2(:,i+1) <= 1.0)

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % step successful
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(iPhase==0)
                    newOrder = 1;
                    iPhase=0;
                    iNoRej=0;
                else
                    iPhase=0;
                    newOrder = 3;
                end

                nef=0;
                ncf=0;
                nsf=0;

                % interpolation between meshpoints
                if(TIME(I+1)==t(i+1))
                    y_interpol(:,I+1)=y(:,i+1);
                    I=I+1;
                elseif(TIME(I+1)<t(i+1))
                    while(TIME(I+1)<t(i+1))
                        if(i>5)
                            nn=5;
                        else
                            nn=i;
                        end
                        y_interpol(:,I+1)=interpol(t,y,TIME(I+1),nn);
                        I=I+1;
                        if(I-1==N)
                            return
                        end
                    end
                end


%            else
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % reject step, step unsucessful
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 iPhase=1;
% 
%                 %restore y,t,y_est,err1
%                 y=y(:,1:i);
%                 t=t(1:i);
%                 y_est=y_est(:,1:i);
%                 err2=err2(:,1:i);
% 
%                 i=i-1;
%                 iNoRej=iNoRej+1;
%                 
%                 if(nef>60)
%                     error('repeated error test failure');
%                 end
% 
%                 %test whether failure was due to corrector iteration  or error test failure
%                 if(converged==0)   % Newton iteration converged
%                     nef = nef+1;
%                     if(nef > 1)
%                         if( nef > 2)
%                             newOrder=2;
%                         else % second error test failure
%                             newOrder=3;
%                         end
%                     else  % first error test failure
%                         newOrder=3;
%                     end
% 
%                 else % Newton iteration failed to converge
%                     if(IER == 0)
%                         ncf=ncf+1;
%                     else
%                         nsf=nsf+1;
%                     end
%                     newOrder=3;
%                 end
% 
%              end
        else
            % interpolation between meshpoints
            if(TIME(I+1)==t(i+1))
                y_interpol(:,I+1)=y(:,i+1);
                I=I+1;
            elseif(TIME(I+1)<t(i+1))
                y_interpol(:,I+1)=interpol(t,y,TIME(I+1),i);
                I=I+1;
            end
        end
    end

    i=i+1;
    
end

