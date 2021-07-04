%% Example Title
% Summary of example objective
tau = 4;
ep = 5e-1; 
h = 0.01; 
t = 2;

s = 0:h:tau;
x = [];

for i = 1:length(s)
    if abs(t-ep-s(i))<ep
        x = [x 1-1/ep*abs(t-ep-s(i))];
    else
        x = [x 0];
    end
end

plot(s,x,'r-','LineWidth',3)
grid on

print -depsc Fig1.eps
  ! epstopdf Fig1.eps
  ! rm Fig1.eps

%% Section 2 Title
% Description of second code block
