% chech order of method

function [C,Ct,ES,XI,eta] = checkOrder(M,c,s,k, order)


if( order ==1 )     % first order


    % first order
    A = M(1:s,1:s);
    U = M(1:s,s+1:s+k);
    B = M(s+1:s+s,1:s );
    V = M(s+1:s+s,s+1:s+k );

    p_max = 3;p=3;
    [dec1,dec2,rho,gamma,gammaD,S,xi,xiD,phi] = algo28(p_max,p,k,s,A,U,B,V);


    for t=0:p_max
        for i=1:s
            if( t==0)
                C(i,t+1)=c(i)^0;
            else
                C(i,t+1)=c(i)^rho(t+1)/gamma(t+1);
            end
        end
        if(t==0)
            E(t+1)=0;
        else
            E(t+1)=1/gamma(t+1);
        end
    end

    ES=comp(E,S,4);

    for i=0:4
        for t=0:4
            if(i==0 && t==0)
                D(i+1,t+1)=1;
            elseif rho(t+1)==i
                D(i+1,t+1)=faculty(rho(t+1))/gamma(t+1);
            else
                D(i+1,t+1)=0;
            end
        end
    end


    for t=0:4
        if( t==0)
            eta(:,t+1)=U*S(:,t+1);
        elseif t>=1 && t<=3
            eta(:,t+1)=A*( eta(:,t) )+U*S(:,t+1);
        elseif t==4
            eta(:,t+1)=A*( eta(:,t-2).*eta(:,t-2) )+U*S(:,t+1);
        end
    end


    X = B*(eta*D);
    Xn(:,1)=[0;0;0];
    Xn(:,2:6)=X;
    Xn=Xn(1:3,1:5);
    XI = Xn+V*S;

elseif( order ==2)  % second order
    % second order
    At = M(1:s,1:s);
    Ut = M(1:s,s+1:s+k);
    A  = M(s+1:s+s,1:s);
    U  = M(s+1:s+s,s+1:s+k);
    Bt = M(2*s+1:2*s+s,1:s );
    Vt = M(2*s+1:2*s+s,s+1:s+k );


    p_max = 3;p=3;
    [dec1,dec2,rho,gamma,gammaD,S,xi,xiD,phi] = algo28(p_max,p,k,s,At,Ut,Bt,Vt,order);


    for t=0:p_max
        for i=1:s
            if( t==0)
                C(i,t+1)=c(i)^0;
                Ct(i,t+1)=c(i)^0;
            else
                C(i,t+1)=c(i)^rho(t+1)/gamma(t+1);
                Ct(i,t+1)=c(i)^rho(t+1)*rho(t+1)/gamma(t+1);
            end
        end
        if(t==0)
            E(t+1)=0;
        else
            E(t+1)=1/gamma(t+1);
        end
    end

    ES=comp(E,S,4);
    XI=0;
    eta=0;
    %
    %     for i=0:4
    %         for t=0:4
    %             if(i==0 && t==0)
    %                 D(i+1,t+1)=1;
    %             elseif rho(t+1)==i
    %                 D(i+1,t+1)=faculty(rho(t+1))/gamma(t+1);
    %             else
    %                 D(i+1,t+1)=0;
    %             end
    %         end
    %     end
    %
    %
    for t=0:4
        if( t==0)
            eta(:,t+1)=Ut*S(:,t+1);
            etaD(:,t+1)=U*S(:,t+1);
        elseif t>=1 && t<=3
            eta(:,t+1)=At*( etaD(:,t)*eta(:,t) )+Ut*S(:,t+1);
            etaD(:,t+1)=A*( etaD(:,t)*eta(:,t) )+U*S(:,t+1);
        elseif t==4
            eta(:,t+1)=At*( eta(:,t-2).*eta(:,t-2) )+Ut*S(:,t+1);
            etaD(:,t+1)=A*( etaD(:,t-2).*etaD(:,t-2) )+U*S(:,t+1);
        end
    end
    %
    %
    %     X = B*(eta*D);
    %     Xn(:,1)=[0;0;0];
    %     Xn(:,2:6)=X;
    %     Xn=Xn(1:3,1:5);
    %     XI = Xn+V*S;

end


