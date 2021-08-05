function [ES,AS,ASd] = qz_transform(E,A,Ad)

[n,m] = size(E);
[~,m1] = size(Ad);

k = floor(m1/n);

if n~=m
    error('Not square matrices')
end

[EE,AA,Q,Z] = qz(E,A,'real');
[ES,AS,QS,ZS] = ordqz(EE,AA,Q,Z,'lhp');

ASd = QS * Ad;
for i = 0:(k-1)
    ASd(:,i*n+1:(i+1)*n) = ASd(:,i*n+1:(i+1)*n) * ZS ;
end