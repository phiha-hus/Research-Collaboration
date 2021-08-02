function [EES,AAS,AASd] = qz_transform(E,A,Ad)

[n,m] = size(E);

if n~=m
    error('Not square matrices')
end

[EE,AA,Q,Z] = qz(E,A);
[EES,AAS,QS,ZS] = ordqz(EE,AA,Q,Z,'lhp');

AASd = QS * Ad * ZS ;
%AASd1 = QS * Ad1 * ZS
%AASd2 = QS * Ad2 * ZS
