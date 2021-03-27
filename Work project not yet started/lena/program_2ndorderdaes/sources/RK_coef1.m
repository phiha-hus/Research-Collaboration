%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge-Kutta coefficients
%
% Returns the Runge-Kutta coefficients
%
%           c| A 
%          --|---
%            | b 
%
% CALL  :
% [A,b,c] = RK_coef(s,i,type)
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
function [A,b,c] = RK_coef1(s,type)


if(strcmp(type,'Radau'))
    if( s==2)
        
        A=[5/12  -1/12
            3/4   1/4];  
        b=[3/4 1/4];
        c=[1/3 1];
    end
    
    if( s==1)
        A=[ 1 ];  
        b=[ 1 ];
        c=[ 1 ];
    end
    
    
elseif(strcmp(type,'Gauss'))  
    if( s==2)
        
        A=[1/4  (3-2*sqrt(3))/12
            (3+2*sqrt(3))/12 1/4];
        b=[1/2  1/2];
        c=[(3-sqrt(3))/6  (3+sqrt(3))/6];
        
    elseif( s==1)
        A=[1/2];
        b=[1];
        c=[1/2]
    else
        error('wrong stage');
    end
    
else
    error('wrong type');
end


