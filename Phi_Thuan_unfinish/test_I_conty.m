% This tests the I-controllability of some first order reformulation
% by minimal extension
%

% Size of the orginal system. May not need to be a square system
d1 = 3; d2 = 3;

M1 = [1 0 0]; 
D2 = [];          D1 = [0 0 0];     D4 = [0 1 0];
K3 = [0 1 0];     K1 = [0 0 0];     K2 = [];
                  K4 = [0 0 0];     K5 = [];


[r2,~] = size(M1);
[r1,~] = size(D2);
[r0,~] = size(K3);
[p1,~] = size(Si1);
[p0,~] = size(Si0);
v = d1 - (r2+r1+r0+p1+p0);


% Si1 = 1; Si0 = [];
B14 = zeros(2*r2+r1+r0,d2);
B5 = [0 1 0];  % coresp. to phi1
B6 = [];        % coresp. to phi0
B7 = zeros(v,d2);

E1 = [eye(r2); zeros(d1,r2)];
E2 = [D1; M1; D2; zeros(size(K3)) ; D4; zeros(p0+v,d2)];
E = [E1 E2]

A1 = [zeros(r2,r2); -eye(r2); zeros(r1+r0+p1+p0+v,r2)];
size(A1);
A2 = [K1; zeros(r2,d2); K2; K3;  K4; K5; zeros(v,d2)];
size(A2);
A = [A1 A2]

B = [B1; B2; B3; B4; B5; B6; B7 ]

M = [E zeros(size(E)) zeros(size(B));A E B]
rank(M) - rank(E) - d2
