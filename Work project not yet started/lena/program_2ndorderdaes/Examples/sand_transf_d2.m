%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = sand_transf_d2(t,y,dy,ddy,type)

if(strcmp(type,'bdf'))   

f= [0  0 0
    0  0 0
    0  0 0];
else

 f= [0  0 0
     0  0 0
     0  0 0];
end