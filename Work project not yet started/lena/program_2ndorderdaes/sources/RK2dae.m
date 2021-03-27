%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge-Kutta method
%
% Solution of second order differential-algebraic equations ( equations of motion)
% of the form
% M q''(t)= f(t,q,q')-G^T(q,t)lambda
%       0 = g(q,t)
%
% with 2-stage Runge-Kutta method.
%
% CALL  :
% [t,h,y,y_2]=RK2dae(funct,y0,y1,t0,a,N,s,RTOL,ATOL,H0,type)
%
% INPUT : funct - name of routine to determin
%                 (Text-String)
%         y0,y1 - initial values
%         t0    - initial time     t\in[t0,t0+a]
%         a     - Intervallaenge   t\in[t0,t0+a]
%         N     - number of steps
%         s     - number of stages, s=2 !
%         RTOL  - relative error tolerance
%         ATOL  - absolute error tolerance
%         H0    - initial stepsize, if H0=0 then the code determines the
%                 initial stepsize
%         type  - type of Runge-Kutta method
%                 type = 'Gauss'
%                        'Radau'
%                        'Lobatto'
%
% OUTPUT: t    - time steps
%         h    - vector of stepsizes
%         y    - approximate solution at times t(i)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,h,y,y_interpol,TIME]=RK2dae(funct,funct_d1,funct_d2,funct_d3,y0,y1,t0,a,N,s,RTOL,ATOL,H0,type,step)

t(1)=t0;

TIME=t0:a/N:t0+a;
TOUT = t0+a;

% initial stepsize
if(H0==0)
    H0 = sign(TOUT-t0)*min(10^(-3)*abs(TOUT-t0),0.5*myNorm(y1,RTOL,ATOL)^(-1));
end

TOUT=t0+a;

h(1)=H0;

[n,m]=size(y0);
[n1,m1]=size(y1);

y(1:n,1)=y0;
y(n+1:n+n1,1)=y1;
y_2(:,1)=y(:,1);
y_1(:,1)=y(:,1);
y_interpol(:,1)=y(:,1);


iNoRej = 0;nef=0;nsf=0;iPhase=0;%newOrder = 0;

i=1;I=1;

%facmax = 1.5;facmin = 0.2;fac = 0.8;
%p=2;

% Integration
while(i<=2000&& t(i)<=TOUT )
    i
    t(i+1) = t(i)+h(i);
    h(i+1) = h(i);

    [y(:,i+1),converged] = RK_step(funct,funct_d1,funct_d2,funct_d3,n,y(:,1:i),i,h(i),t(1:i),s,type,RTOL,ATOL);

    % Error Estimate - Extrapolationsmethode
    if(i<=4)
        y_est(:,i+1)=extrapol(t,y,i-1);
        err2(:,i+1)=myNorm(y_est(:,i+1)-y(:,i+1),RTOL,ATOL);
    end
    if(i>4)
        y_est(:,i+1)=extrapol(t,y,3);
        err2(:,i+1)=myNorm(y_est(:,i+1)-y(:,i+1),RTOL,ATOL);
    end

    hmin = 4*eps*max(abs(t(i+1)),abs(TOUT));

    EST=err2(:,i+1);

    if( step ==1)

        if( i>5)
            if( err2(:,i+1) <= 1.0)

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % step successful
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(iPhase==0)
                    h(i+1)=2*h(i);
                    iPhase=0;
                    iNoRej=0;
                else
                    r=(2*EST+0.0001)^(-1/(1+1));
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
                end

                nef=0;
                %ncf=0;
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
                            r=0.25;
                            h(i+1)=r*h(i+1);
                            if(abs(h(i+1))<hmin)
                                error('minimal stepsize');
                            end
                        else % second error test failure
                            r=0.25;
                            h(i+1)=r*h(i+1);
                            if(abs(h(i+1))<hmin)
                                error('minimal stepsize');
                            end
                        end
                    else  % first error test failure
                        r=0.9*(2.0*EST+0.0001)^(-1/(1+1));
                        r=max(0.25,min(0.9,r));
                        h(i+1)=r*h(i+1);
                        if(abs(h(i+1))<hmin)
                            error('minimal stepsize');
                        end
                    end

                else % Newton iteration failed to converge
                    nsf=nsf+1;
                    r=0.25;
                    h(i+1)=r*h(i+1);
                    if(abs(h(i+1))<hmin || nsf >3 )
                        error('minimal stepsize');
                    end

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
        h(i+1)=h(i);
        if( i>10)
   %         if( err2(:,i+1) <= 1.0)

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % step successful
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(iPhase==0)
                    iPhase=0;
                    iNoRej=0;
                else
                    iPhase=0;
                end

                nef=0;
                %ncf=0;
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

%             else
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
%                 %test whether failure was due to corrector iteration  or error test failure
%                 if(converged==0)   % Newton iteration converged
%                     nef = nef+1;
%                 else % Newton iteration failed to converge
%                     nsf=nsf+1;
%                 end
%             end
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

