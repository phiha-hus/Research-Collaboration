% runme 
clear;
figure(1);clf;
            
a = 5;
N = 10; 
h = 0.001;
            
RTOL=1e-4*ones(15,1);
ATOL=1e-4*ones(15,1);
            
step = 2; % constant stepsize
        
[t0,y0]=feval('example1_para');
n=length(y0);
[t1,h1,y1,y1_est,y1_i]=bdfk1_dae('example1','example1_d1','example1_d2',y0,t0,a,N,RTOL,ATOL,h,step,3); 

% [t0,y0]=feval('example1_para');
% n=length(y0);
% [t2,h2,y2,y2_est,y2_i]=bdfk1_dae('example2','example2_d1','example2_d2',y0,t0,a,N,RTOL,ATOL,h,step,3); 

figure(1);
subplot(1,3,1)
plot(t1,y1(1,:),'b');hold on;
%plot(t2,y2(1,:),'r');hold on;
xlabel('t');
ylabel('x');
axis tight;
subplot(1,3,2)
plot(t1,y1(2,:),'b');
%plot(t2,y2(2,:),'r');
xlabel('t');
ylabel('y');
axis tight;
subplot(1,3,3)
plot(t1,y1(3,:),'b');
%plot(t2,y2(3,:),'r');
axis tight;
xlabel('t');
ylabel('L');        
    