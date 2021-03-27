%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exact solution of example sand
%
% CALL  : ye = example1_exact(t)
%
%
% INPUT : t   - current time t\in[t0,t0+a]
%
% OUTPUT: ye   - exact solution at t
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ye = sand_exact(t)

ye=[sin(t^2)
    cos(t^2)
    -4*t^2
    2*t*cos(t^2)
   -2*t*sin(t^2)];
    %-8*t];


