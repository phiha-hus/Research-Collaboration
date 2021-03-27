function [trig_flag,t_trig,x_trig,idx] = find_trig_point(T,x,sigma)
% Given the time T and the corresponding discrete signal x
% We want to find the first triggering point that match the condition
% norm(x-x(:,1))./norm(x) >= sigma
% idx : index of the first point t_{k+1} where the condition is matched
% Copyright Phi Ha 18/06/2020

[~,nT] = size(T);
[~,n] = size(x);

if nT ~= n
    error('sizes of time and signal x do not match')
end

abs_err = x - kron(ones(1,n),x(:,1));

rel_err = zeros(1,n);

for i = 1:n
    rel_err(i) = norm(abs_err(:,i))/norm(x(:,i));
end

ID = find( rel_err >= sigma) ; 
size(ID)

if size(ID) ~= 0
    trig_flag = 1;
    idx = ID(1);
    t_trig = T(idx);
    x_trig = x(:,idx);
else
    trig_flag = 0;
    idx = [];
    t_trig = [];
    x_trig = [];
end
    
        

        



