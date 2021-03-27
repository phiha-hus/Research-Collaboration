%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation of Predictorpolynomial
%
%
% called by phi_star
%         
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = psi(i,N,H) %N=n+1

sum=0;
for c=0:i-1
    sum=sum+H(N-c);
end
p = sum ;

