% Fanbin Bu and Yimin Wei, The algorithm for computing the ...
% ... Drazin inverses of two-variable polynomial matrices, 
% Applied mathematics and computation 147.3 (2004): 805-836.
%
function [Ad] = Drazin_inverse(a)
%-----------------------------------------
% Compute the Drazin Inverse of a matrix 'a' using the limited 
% algorithm.
% Need computing the index of 'a'.
[m,n] = size(a);
if m~= n
    disp('Matrix is must be square!')
end
%-----------------------------------------
% Computer the index of A and note r = rank(A^k).
[k,r,ak] = index(a);

% Initialize F1 and g1
F = eye(n);
g = -trace(a^(k+1));

% Begin the loop
ak1 = ak * a;
for i = 1:r-1
    F = ak1*F + g*eye(n);
    g = -trace(ak1*F)/(i+1);    
end

Ad = -F/g * ak;
end

function [k,r,ak] = index(a)
%To compute the index of 'a'.
k = 0;
n = size(a);
ak = a^2 ;
r = rank(a);
rk = rank(ak) ;

if rank(a) == n
    disp('Warning: non-singular matrix')
    a
end

while rk ~= r
    a = ak; r = rk;
    ak = ak * a; rk = rank(ak);
    k = k+1;    
end

end