function ye = sand_transf_exact(t)

ye=[sin(t^2)-4*t^2
    cos(t^2)
    -4*t^2
    2*t*cos(t^2)-8*t
   -2*t*sin(t^2)];
    %-8*t];