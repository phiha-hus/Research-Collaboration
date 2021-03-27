function y=faculty(x)

if( x<0)
    error('wrong argumnet');
end

if(x==0)
    y=1;
else
    y=x*faculty(x-1); 
end