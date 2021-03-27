%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition of matrix for example 1 - GGL
%
% CALL  : M=example1_RK(t,T,h,y,yi,K,a,sel)
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
%                 sel=0 - dimension of [q,q']
%                 sel=1 - evaluation of F
%                 sel=2 - evaluation of DF
%
% OUTPUT: M   - matrix or vector 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=example1_rot(t,y,dy,ddy)

    f= [ddy(1)+ y(1)+2*y(1)*dy(1)*dy(2)+2*y(1)*y(3)
        ddy(2)+ dy(1)-2*y(1)*y(2)^2+2*y(2)*y(3)
        y(1)^2+y(2)^2-1+2*y(1)*dy(1)*dy(2)+2*y(1)*y(3)];