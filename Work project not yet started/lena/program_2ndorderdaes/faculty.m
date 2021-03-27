function y=faculty(x)

if(x<=0)
    if(x==0)
        y=1;
    else
        y=0;
    end
else
    y=x*faculty(x-1); 
end