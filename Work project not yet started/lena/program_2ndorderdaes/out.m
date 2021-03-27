figure(1);
subplot(2,3,1); 
%plot(t,ye(1,:),'r'); hold on; 
plot(t,y_bdf(1,:),'b--'); hold on;
plot(t,y_bdf1(1,:),'b');
xlabel('t'); 
ylabel('x');
axis tight;

subplot(2,3,2); 
%plot(t,ye(2,:),'r'); hold on; 
plot(t,y_bdf(2,:),'b--'); hold on;
plot(t,y_bdf1(2,:),'b');
xlabel('t'); 
ylabel('y');
axis tight;

subplot(2,3,3);
%plot(t,ye(3,:),'r'); hold on; 
plot(t,y_bdf(3,:),'b--'); hold on;
plot(t,y_bdf1(5,:),'b');
xlabel('t'); 
ylabel('\lambda');
axis tight;

% subplot(2,3,4);
% %plot(t,ye(4,:),'r'); hold on; 
% plot(t,y_bdf(4,:),'b--'); hold on;
% plot(t,y_bdf1(3,:),'g');
% xlabel('t'); 
% ylabel('x\prime');
% axis tight;

% subplot(2,3,5);
% %plot(t,ye(5,:),'r'); hold on; 
% plot(t,y_bdf(5,:),'b--'); hold on;
% plot(t,y_bdf1(4,:),'g');
% xlabel('t'); 
% ylabel('y\prime');
% axis tight;


subplot(2,3,4);
semilogy(t,err_bdf(1,:),'b--');hold on;
semilogy(t,err_bdf1(1,:),'b');
xlabel('t');
ylabel('error x');
axis tight;
        
subplot(2,3,5);    
semilogy(t,err_bdf(2,:),'b--');hold on;
semilogy(t,err_bdf1(2,:),'b');
xlabel('t');
ylabel('error y');
axis tight;

subplot(2,3,6);
semilogy(t,err_bdf(3,:),'b--');hold on;
semilogy(t,err_bdf1(5,:),'b');
xlabel('t');
ylabel('error \lambda');
axis tight;
        
% subplot(2,3,9);
% semilogy(t,err_bdf(4,:),'b--');hold on;
% semilogy(t,err_bdf1(3,:),'g');
% xlabel('t');
% ylabel('error x\prime');
% axis tight;
% 
% subplot(2,3,10);
% semilogy(t,err_bdf(5,:),'b--');hold on;
% semilogy(t,err_bdf1(4,:),'g');
% xlabel('t');
% ylabel('error y\prime');
% axis tight;