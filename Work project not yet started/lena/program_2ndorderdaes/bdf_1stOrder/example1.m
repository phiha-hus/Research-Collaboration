function f = example1(t,y,dy)

f= [dy(1)-y(4)
    dy(2)-y(5)
    dy(4)-2*y(2)-y(1)*y(3)
    dy(5)+2*y(1)-y(2)*y(3)
    y(1)^2+y(2)^2-1];
