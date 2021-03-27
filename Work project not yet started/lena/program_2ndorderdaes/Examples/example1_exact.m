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
function ye = example1_exact(t)

ye=[sin(t)
    cos(t)
    sin(t)*cos(t)
    cos(t)
    -sin(t)];
    %cos(t)^2-sin(t)^2];
