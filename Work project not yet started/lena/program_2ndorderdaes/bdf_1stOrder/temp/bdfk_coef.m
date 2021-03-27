%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BDFk-coefficients
%
% returns bdf-coefficients for order k
%
% CALL  : a=bdfk_coef(k,j,l,T)
%
% INPUT :    k		- order of method
%            j,l    - BDF(l,j,k) method
%			 T		- time interval T=[t(i),..,t(i-k-2)]
%            K      - order of last taken steps
%				
%
% OUTPUT: a    - coefficients of the method a=[a1,..,ak];
%         
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [a]=bdfk_coef(K,T)%(k,j,l,T)	% K=(k,j,l,...)	% BDF(l,j,k)
function alpha=bdfk_coef(K,T)%(k,j,l,T)		% BDF(l,j,k)

% k=K(1);
% j=K(2);
% l=K(3);
% 
% if( k==1 )	% method of order 1
%     if( j==1)
%         a = 2/(T(1)-T(3));		
%     elseif(j==2)
%         a = 2/(T(1)-T(2));    
%     else
%         error('wrong arguments')
%     end
%     
% elseif( k==2)	% method of order 2
%     if( j==1 )
%         if( l==1)		% BDF(1,1,2)
%             a=[2*(3*T(1)-T(2)-T(3)-T(4))/((2*T(1)-T(2)-T(3))*(T(1)-T(4))),-2*(2*T(1)-T(2)-T(3))/((T(2)-T(4))*(T(1)-T(4))) ];        
%             
%         elseif( l==2) % BDF(2,1,2)
%             b5=2*(2*T(1)-T(2)-T(3))/((T(1)-T(4))*(T(2)-T(3))-(T(2)-T(5))*(T(3)-T(4)));
%             b4=(b5*(T(2)-T(3))+2)/(2*T(1)-T(2)-T(3));
%             a=[b4,-b5];
%             
%         elseif( l==3)	% BDF(3,1,2)    
%             a2 = 2*(2*T(1)-T(2)-T(3));
%             a=[ 2*(3*T(1)-T(2)-2*T(3))/((2*T(1)-T(2)-T(3))*(T(1)-T(3))), -a2/((T(1)-T(3))*(T(2)-T(3)))];    
%         else
%             error('wrong argument');
%         end
%         
%     elseif( j==2 )
%         if(l==1)		% BDF(1,2,2)
%             b1=2*(3*T(1)-T(2)-T(3)-T(4))/((2*T(1)-T(2)-T(3))*(T(1)-T(4)));
%             b3=2*(b1*(T(1)-T(2))-1)/(2*T(2)-T(3)-T(4));
%             a=[b1,-b3];
%         elseif(l==2)	% BDF(2,2,2)        
%             a2=(T(1)-T(3))*(T(2)-T(3))*(2*(3*T(1)-T(2)-T(3)-T(4))*(T(1)-T(2))-(2*T(1)-T(2)-T(3))*(T(1)-T(4)))/((2*T(1)-T(2)-T(3))*(T(1)-T(4))*(T(2)-T(3))-(T(1)-T(2))*(T(2)-T(5))*(T(3)-T(4)));
%             a=[1/(T(1)-T(2))+a2/((T(1)-T(3))*(T(1)-T(2))) , -a2/((T(1)-T(3))*(T(2)-T(3)))];
%         else
%             error('wrong argument');
%         end
%         
%     else
%         error('wrong argument');
%     end
%     
% else
%     error('wrong arguments')
% end






% Allgemein A*alpha = b
%K=[k,j,l];
k=K(1);

%b = zeros(k,1);
%y_1 =zeros(k,k+1);
%A = zeros(k,k);


for i=1:k
    b(i)=faculty(i+1)/faculty(i-1)*T(1)^(i-1);
end

for i=1:k
    for j=1:k

        y_0(i,:)=T.^(i+1);
        
        for c=1:j+1
            y_1(i,c)=bdf1( y_0(i,c:c+K(c)),T(c:c+K(c)),K(c));
        end

        A(i,j)=divDiff(T(1:j+1),y_1(i,1:j+1),j);
    end
end

alpha = A\b';
alpha = alpha';



