%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM-Verfahren, variable-step variable-order GLM
%
% Loesen der differentiell-algebraischen Gleichung 2. Ordung der Form
% M y''(t)=D y'(t) + S y(t) + f(t,y,y')
%
% mit GLM der Ordnung Ordnung k
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
function [t,h,y,yp]=glm1dae(funct,funct_start,funct_d1,funct_d2,y0,t0,a,N,p,q,k,s,RTOL,ATOL,H0,step)

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

iNoRej = 0;nef=0;nsf=0;iPhase=0;

%facmax = 1.5;
%facmin = 0.2;
%fac = 0.8;

% integration
i=1;I=1;

%[y(:,1)] = starting1(funct,funct_d1,funct_d2,n,p,q,k,s,h,t0,y0,RTOL,ATOL);
%y(:,1)=[1;-1;2;-10*h(1);10*h(1);-20*h(1);100*h(1)^2;-100*h(1)^2;200*h(1)^2];
y(:,1)=feval(funct_start,t0,h(1),s);
yp(:,1)  = y(n+1:2*n,1)/h(1);

Y(:,1)    = kron(ones(s,1),eye(n,n))*y0;
dY(:,1)   = kron(ones(s,1),eye(n,n))*y0;

% y_interpol(:,1) = y(:,1);
y2(1:n,1) = y(1:n,1);
y2(n+1:2*n,1) = yp(:,1);

while(i<700 && t(i)<=TOUT)
    i
    t(i+1) = t(i)+h(i);
    h(i+1) = h(i);

    % GLM coefficients
    [M,c] = glm1coef(p);
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     A = M(1:s,1:s);
%     U = M(1:s,s+1:s+k);
%     B = M(s+1:s+s,1:s );
%     V = M(s+1:s+s,s+1:s+k );
%     stabilityRegion(A,A,U,U,V,B,s);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %internal stages -> Newton
    [Y(:,i+1),dY(:,i+1),converged]= glm1_internal(funct,funct_d1,funct_d2,dY(:,i),y(:,i),n,M,c,h,i,t(i),p,s,k,RTOL,ATOL,0 );

    % output quantities
    [y(:,i+1)] =  glm1_external(y(:,i),dY(:,i+1),n,M,h,i,p,s,k );
    yp(:,i+1)  =  y(n+1:2*n,i+1)/h(i);

    y2(1:n,i+1) = y(1:n,i+1);
    y2(n+1:2*n,i+1) = yp(:,i+1);

    % Error Estimate - Extrapolationsmethode
    if(i<=4)
        y_est(:,i+1)=extrapol(t,y2,i-1);
        err2(:,i+1)=myNorm(y_est(:,i+1)-y2(:,i+1),RTOL,ATOL);
    end
    if(i>4)
        y_est(:,i+1)=extrapol(t,y2,3);
        err2(:,i+1)=myNorm(y_est(:,i+1)-y2(:,i+1),RTOL,ATOL);
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
        if( i>5)
            if( err2(:,i+1) <= 1.0)

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

            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % reject step, step unsucessful
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                iPhase=1;

                restore y,t,y_est,err1;
                y=y(:,1:i);
                t=t(1:i);
                y_est=y_est(:,1:i);
                err2=err2(:,1:i);

                i=i-1;
                iNoRej=iNoRej+1;

                %test whether failure was due to corrector iteration  or error test failure
                if(converged==0)   % Newton iteration converged
                    nef = nef+1;
                else % Newton iteration failed to converge
                    nsf=nsf+1;
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
    end


    i=i+1;
end


















