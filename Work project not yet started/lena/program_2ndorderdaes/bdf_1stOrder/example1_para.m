%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% initial values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,y0,y1] = example1_para()


t=0;
y0=[ sin(t); cos(t);sin(t)*cos(t);cos(t);-sin(t)];
y1=[cos(t);-sin(t);(cos(t))^2-(sin(t))^2;-sin(t);-cos(t)];