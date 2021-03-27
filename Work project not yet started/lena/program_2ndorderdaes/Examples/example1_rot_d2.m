function f=example1_rot_d2(t,y,dy,ddy,type)

  
if(strcmp(type,'bdf'))   
      
f = [2*y(1)*dy(2)  2*y(1)*dy(1) 
      1  0  
     2*y(1)*dy(2)  2*y(1)*dy(1) ];
  
else
    
    f = [2*y(1)*dy(2)  2*y(1)*dy(1) 0
      1  0  0
      2*y(1)*dy(2)  2*y(1)*dy(1) 0];
  
end