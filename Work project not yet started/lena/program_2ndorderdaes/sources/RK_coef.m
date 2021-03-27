%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge-Kutta coefficients
%
% Returns the Runge-Kutta coefficients
%
%           c| A |At
%          --|---|---
%            | b |bt
%
% CALL  :
% [A,At,b,bt,c] = RK_coef(s,i,type)
%
% INPUT : s     - number of stages, s=2 !
%         i     - number of current step
%         type  - type of Runge-Kutta method
%                 type = 'Gauss'
%                        'Radau'
%              
%
% OUTPUT: A    - A= [a11 a12
%                    a21 a22]
%         At   - At=[at11 at12
%                    at21 at22]
%         b    - b= [b1 b2]
%         bt   - bt=[bt1 bt2]
%         c    - c= [c1 c2]'
%         
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,At,b,bt,c] = RK_coef(s,type)


if(strcmp(type,'Radau'))
    if( s==2)
        
        A=[5/12  -1/12
            3/4   1/4]; 
        At=[1/9  -1/18
            1/2   0  ]; 
        b=[3/4 1/4];
        bt=[1/2 0];
        c=[1/3 1];
        
    elseif(s==1)
        A=[1];
        b=[1];
        c=[1];
        At=[1];
        bt=[1];
        
    elseif( s==3)
        % Radau5
        A=[(88-7*sqrt(6))/360  (296-169*sqrt(6))/1800 (-2+3*sqrt(6))/225
           (296+169*sqrt(6))/1800  (88+7*sqrt(6))/360 (-2-3*sqrt(6))/225
           (16-sqrt(6))/36 (16+sqrt(6))/36  1/9]; 
        At=[0.0218   -0.0199    0.0100
            0.1772    0.0382   -0.0074
            0.3180    0.1820         0 ]; 
        b=[(16-sqrt(6))/36 (16+sqrt(6))/36  1/9];
        bt=[0.3180    0.1820         0];
        c=[(4-sqrt(6))/10 (4+sqrt(6))/10 1];
        
    else
        error('wrong stage');
    end
    
    
elseif(strcmp(type,'Gauss'))  
    if( s==2)
        
        A=[1/4  (3-2*sqrt(3))/12
            (3+2*sqrt(3))/12 1/4];
        At=[1/24  (3-2*sqrt(3))/24
            (3+2*sqrt(3))/24  1/24];
        b=[1/2  1/2];
        bt=[(3+sqrt(3))/12  (3-sqrt(3))/12];
        c=[(3-sqrt(3))/6  (3+sqrt(3))/6];
        
    elseif( s==1)
        A = [1/2];
        At = [1/4];
        b=[1];
        bt=[1/2];
        c=[1/2];
    else
        error('wrong stage');
    end
    
else
    error('wrong type');
end


