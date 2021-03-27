function ye = sand_1_start(t0,h,s)


if( s==3)
ye=[sin(t0^2)
    cos(t0^2)
    -4*t0^2
    2*t0*cos(t0^2)
    -2*t0*sin(t0^2)
    h*2*t0*cos(t0^2)
    -h*2*t0*sin(t0^2)
    -h*8*t0
    h*(2*cos(t0^2)-4*t0^2*sin(t0^2))
    -h*(2*sin(t0^2)+4*t0^2*cos(t0^2))
    h^2*(2*cos(t0^2)-4*t0^2*sin(t0^2))
    -h^2*(2*sin(t0^2)+4*t0^2*cos(t0^2))
    -h^2*8
    -h^2*(4*t0*sin(t0^2)+8*t0*sin(t0^2)+8*t0^3*cos(t0^2))
    -h^2*(4*t0*cos(t0^2)+8*t0*cos(t0^2)-8*t0^3*sin(t0^2))];
end