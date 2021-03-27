function f=example1_rot_d1(t,y,dy,ddy)

f = [ 1+2*dy(1)*dy(2)+2*y(3)  0  2*y(1)
      -2*y(2)^2  -4*y(1)*y(2)+2*y(3)  2*y(2)
      2*y(1)+2*dy(1)*dy(2)+2*y(3)  2*y(2)  2*y(1)];