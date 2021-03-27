%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exact solution of example 1
%
% CALL  : ye = example1_exact(t)
%
%
% INPUT : t   - current time t\in[t0,t0+a]
%
% OUTPUT: ye   - exact solution at t
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ye = example1_start(t0,h,s)

if( s==3)
ye=[sin(t0)
    cos(t0)
    sin(t0)*cos(t0)
    h*cos(t0)
    -h*sin(t0)
    h*(cos(t0)^2-sin(t0)^2)
    -h^2*sin(t0)
    -h^2*cos(t0)
    -h^2*4*sin(t0)*cos(t0)];
end