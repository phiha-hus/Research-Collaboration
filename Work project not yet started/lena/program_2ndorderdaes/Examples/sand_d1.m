%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition of matrix for example [sand]
% 
% CALL  : M=example1(t,T,h,y,yi,K,a,sel)
% 
% 
% INPUT : t   - current time t\in[t0,t0+a]
%         T   - time steps T=[t0,t1,..,t]
%         h   - current stepsize
%         y   - vector of current values
%         yi  - vector of previous values
%         K   - current order
%         a   - vector of bdf coefficients
%         sel - selection
%                 sel=1 - dimension of [q,r,w,s,lambda]
%                 sel=2 - evaluation of F
%                 sel=3 - evaluation of DF
% 
% OUTPUT: M   - matrix or vector 
% 
% Author: Lena Wunderlich
% Date:   27.10.2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = sand_d1(t,y,dy,ddy)

f= [-y(3)  -2  -y(1)
    2     -y(3) -y(2)
    2*y(1) 2*y(2) 0];

