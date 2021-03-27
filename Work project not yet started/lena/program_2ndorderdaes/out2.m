figure(1);
subplot(2,3,1); 
%plot(t,ye(1,:),'r'); hold on; 
plot(tk,y_rk(1,:),'b--'); hold on;
plot(tk1,y_rk1(1,:),'b');
xlabel('t'); 
ylabel('x');
axis tight;

subplot(2,3,2); 
%plot(t,ye(2,:),'r'); hold on; 
plot(tk,y_rk(2,:),'b--'); hold on;
plot(tk1,y_rk1(2,:),'b');
xlabel('t'); 
ylabel('y');
axis tight;

subplot(2,3,3);
%plot(t,ye(3,:),'r'); hold on; 
plot(tk,y_rk(3,:),'b--'); hold on;
plot(tk1,y_rk1(5,:),'b');
xlabel('t'); 
ylabel('\lambda');
axis tight;

% subplot(2,5,4);
% plot(t,ye(4,:),'r'); hold on; 
% plot(tk,y_rk(4,:),'b--'); hold on;
% plot(tk1,y_rk1(3,:),'g');
% xlabel('t'); 
% ylabel('x\prime');
% axis tight;
% 
% subplot(2,5,5);
% plot(t,ye(5,:),'r'); hold on; 
% plot(tk,y_rk(5,:),'b--'); hold on;
% plot(tk1,y_rk1(4,:),'g');
% xlabel('t'); 
% ylabel('y\prime');
% axis tight;


subplot(2,3,4);
semilogy(tk,err_rk(1,:),'b--');hold on;
semilogy(tk1,err_rk1(1,:),'b');
xlabel('t');
ylabel('error x');
axis tight;
        
subplot(2,3,5);    
semilogy(tk,err_rk(2,:),'b--');hold on;
semilogy(tk1,err_rk1(2,:),'b');
xlabel('t');
ylabel('error y');
axis tight;

subplot(2,3,6);
semilogy(tk,err_rk(3,:),'b--');hold on;
semilogy(tk1,err_rk1(5,:),'b');
xlabel('t');
ylabel('error \lambda');
axis tight;
        
% subplot(2,5,9);
% semilogy(tk,err_rk(4,:),'b--');hold on;
% semilogy(tk1,err_rk1(3,:),'g');
% xlabel('t');
% ylabel('error x\prime');
% axis tight;
% 
% subplot(2,5,10);
% semilogy(tk,err_rk(5,:),'b--');hold on;
% semilogy(tk1,err_rk1(4,:),'g');
% xlabel('t');
% ylabel('error y\prime');
% axis tight;