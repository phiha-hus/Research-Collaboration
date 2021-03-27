%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = sand_transf_d3(t,y,dy,ddy,type)


if(strcmp(type,'bdf'))
    
f= [1  0 -1
    0  1  0
    0  0  0];

else

 f= [1  0 -1
     0  1  0
     0  0  0];
end
