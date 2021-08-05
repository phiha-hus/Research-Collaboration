clear all; close all; clc

A = [-2 1;1 -1]
B = [1;1]

setlmis([])

P = lmivar(1,[2,1]);
Q = lmivar(1,[2,1]);
R = lmivar(1,[1,1]);

lmiterm([-1 1 1 P],1,1)

lmiterm([-2 1 1 Q],1,1)

lmiterm([-3 1 1 R],1,1)

lmiterm([4 1 1 P],A',1,'s');
lmiterm([4 1 1 Q],1,1);


lmiterm([4 1 2 P],1,B);

lmiterm([4 2 2 R],1,-1);

LMISYS = getlmis;

[T_min,X_feas] = feasp(LMISYS);

Q = dec2mat(LMISYS,X_feas,Q)
P = dec2mat(LMISYS,X_feas,P)
R = dec2mat(LMISYS,X_feas,R)


K = inv(R) * B' * P