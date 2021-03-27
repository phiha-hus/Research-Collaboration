%% OUTPUT

% figure(1);
% subplot(2,3,1); 
% %plot(t,ye(1,:),'r'); hold on; 
% plot(t_bdf2,y_bdf2(1,:),'b--'); hold on;
% plot(t_bdf1,y_bdf1(1,:),'b');
% xlabel('t'); 
% ylabel('x');
% axis tight;
% 
% subplot(2,3,2); 
% %plot(t,ye(2,:),'r'); hold on; 
% plot(t_bdf2,y_bdf2(2,:),'b--'); hold on;
% plot(t_bdf1,y_bdf1(2,:),'b');
% xlabel('t'); 
% ylabel('y');
% axis tight;
% 
% subplot(2,3,3);
% %plot(t,ye(3,:),'r'); hold on; 
% plot(t_bdf2,y_bdf2(3,:),'b--'); hold on;
% plot(t_bdf1,y_bdf1(3,:),'b');
% xlabel('t'); 
% ylabel('\lambda');
% axis tight;
% 
% 
% subplot(2,3,4);
% semilogy(t_bdf2,err_bdf2(1,:),'b--');hold on;
% semilogy(t_bdf1,err_bdf1(1,:),'b');
% xlabel('t');
% ylabel('error x');
% axis tight;
%         
% subplot(2,3,5);    
% semilogy(t_bdf2,err_bdf2(2,:),'b--');hold on;
% semilogy(t_bdf1,err_bdf1(2,:),'b');
% xlabel('t');
% ylabel('error y');
% axis tight;
% 
% subplot(2,3,6);
% semilogy(t_bdf2,err_bdf2(3,:),'b--');hold on;
% semilogy(t_bdf1,err_bdf1(3,:),'b');
% xlabel('t');
% ylabel('error \lambda');
% axis tight;
        


% figure(2);
% subplot(2,3,1); 
% %plot(t,ye(1,:),'r'); hold on; 
% plot(t_rk2,y_rk2(1,:),'b--'); hold on;
% plot(t_rk1,y_rk1(1,:),'b');
% xlabel('t'); 
% ylabel('x');
% axis tight;
% 
% subplot(2,3,2); 
% %plot(t,ye(2,:),'r'); hold on; 
% plot(t_rk2,y_rk2(2,:),'b--'); hold on;
% plot(t_rk1,y_rk1(2,:),'b');
% xlabel('t'); 
% ylabel('y');
% axis tight;
% 
% subplot(2,3,3);
% %plot(t,ye(3,:),'r'); hold on; 
% plot(t_rk2,y_rk2(3,:),'b--'); hold on;
% plot(t_rk1,y_rk1(3,:),'b');
% xlabel('t'); 
% ylabel('\lambda');
% axis tight;
% 
% 
% subplot(2,3,4);
% semilogy(t_rk2,err_rk2(1,:),'b--');hold on;
% semilogy(t_rk2,err_rk1(1,:),'b');
% xlabel('t');
% ylabel('error x');
% axis tight;
%         
% subplot(2,3,5);    
% semilogy(t_rk2,err_rk2(2,:),'b--');hold on;
% semilogy(t_rk1,err_rk1(2,:),'b');
% xlabel('t');
% ylabel('error y');
% axis tight;
% 
% subplot(2,3,6);
% semilogy(t_rk2,err_rk2(3,:),'b--');hold on;
% semilogy(t_rk1,err_rk1(3,:),'b');
% xlabel('t');
% ylabel('error \lambda');
% axis tight;


% figure(3);
% subplot(2,3,1); 
% %plot(t,ye(1,:),'r'); hold on; 
% plot(t_glm2,y_glm2(1,:),'b--'); hold on;
% plot(t_glm1,y_glm1(1,:),'b');
% xlabel('t'); 
% ylabel('x');
% axis tight;
% 
% subplot(2,3,2); 
% %plot(t,ye(2,:),'r'); hold on; 
% plot(t_glm2,y_glm2(2,:),'b--'); hold on;
% plot(t_glm1,y_glm1(2,:),'b');
% xlabel('t'); 
% ylabel('y');
% axis tight;
% 
% subplot(2,3,3);
% %plot(t,ye(3,:),'r'); hold on; 
% plot(t_glm2,y_glm2(3,:),'b--'); hold on;
% plot(t_glm1,y_glm1(3,:),'b');
% xlabel('t'); 
% ylabel('\lambda');
% axis tight;
% 
% 
% subplot(2,3,4);
% semilogy(t_glm2,err_glm2(1,:),'b--');hold on;
% semilogy(t_glm1,err_glm1(1,:),'b');
% xlabel('t');
% ylabel('error x');
% axis tight;
%         
% subplot(2,3,5);    
% semilogy(t_glm2,err_glm2(2,:),'b--');hold on;
% semilogy(t_glm1,err_glm1(2,:),'b');
% xlabel('t');
% ylabel('error y');
% axis tight;
% 
% subplot(2,3,6);
% semilogy(t_glm2,err_glm2(3,:),'b--');hold on;
% semilogy(t_glm1,err_glm1(3,:),'b');
% xlabel('t');
% ylabel('error \lambda');
% axis tight;


% % Circuit
    figure(4);
        subplot(2,3,1)
        plot(t,ye(1,:),'r'); hold on; 
        plot(t_bdf2,y_bdf2(1,:),'b'); 
        plot(t_rk2,y_rk2(1,:),'b--'); 
        plot(t_glm2,y_glm2(1,:),'b:'); 
   
        xlabel('t'); 
        ylabel('v_L');
        axis tight;
        subplot(2,3,2)
        plot(t,ye(2,:),'r'); hold on; 
        plot(t_bdf2,y_bdf2(2,:),'b'); 
        plot(t_rk2,y_rk2(2,:),'b--'); 
        plot(t_glm2,y_glm2(2,:),'b:'); 
        xlabel('t');
        ylabel('v_C');
        axis tight;  
        subplot(2,3,3)
        plot(t,ye(3,:),'r'); hold on; 
        plot(t_bdf2,y_bdf2(3,:),'b'); 
        plot(t_rk2,y_rk2(3,:),'b--'); 
        plot(t_glm2,y_glm2(3,:),'b:'); 
        xlabel('t');
        ylabel('v_R');
        axis tight;  
  
        
        subplot(2,3,4)
        semilogy(t_bdf2,err_bdf2(1,:),'b'); hold on; 
        semilogy(t_rk2,err_rk2(1,:),'b--');
        semilogy(t_glm2,err_glm2(1,:),'b:');
        xlabel('t');
        ylabel('error in v_L');
        axis tight;
        subplot(2,3,5)
        semilogy(t_bdf2,err_bdf2(2,:),'b'); hold on; 
        semilogy(t_rk2,err_rk2(2,:),'b--');
        semilogy(t_glm2,err_glm2(2,:),'b:');
        xlabel('t');
        ylabel('error in v_C');
        axis tight;
        subplot(2,3,6)
        semilogy(t_bdf2,err_bdf2(3,:),'b'); hold on; 
        semilogy(t_rk2,err_rk2(3,:),'b--');
        semilogy(t_glm2,err_glm2(3,:),'b:');
        xlabel('t');
        ylabel('error in v_R');
        axis tight;
        
%            figure(5);
%         subplot(2,3,1)
%         %plot(t,ye(1,:),'r'); hold on; 
%         plot(t_bdf1,y_bdf1(1,:),'b');  hold on; 
%         plot(t_rk1,y_rk1(1,:),'b--'); 
%         plot(t_glm1,y_glm1(1,:),'b:'); 
%    
%         xlabel('t'); 
%         ylabel('v_L');
%         axis tight;
%         subplot(2,3,2)
%         %plot(t,ye(2,:),'r'); hold on; 
%         plot(t_bdf1,y_bdf1(2,:),'b');  hold on; 
%         plot(t_rk1,y_rk1(2,:),'b--'); 
%         plot(t_glm1,y_glm1(2,:),'b:'); 
%         xlabel('t');
%         ylabel('v_C');
%         axis tight;  
%         subplot(2,3,3)
%         %plot(t,ye(3,:),'r'); hold on; 
%         plot(t_bdf1,y_bdf1(3,:),'b');  hold on; 
%         plot(t_rk1,y_rk1(3,:),'b--'); 
%         plot(t_glm1,y_glm1(3,:),'b:'); 
%         xlabel('t');
%         ylabel('v_R');
%         axis tight;  
%   
%         
%         subplot(2,3,4)
%         semilogy(t_bdf1,err_bdf1(1,:),'b'); hold on; 
%         semilogy(t_rk1,err_rk1(1,:),'b--');
%         semilogy(t_glm1,err_glm1(1,:),'b:');
%         xlabel('t');
%         ylabel('error in v_L');
%         axis tight;
%         subplot(2,3,5)
%         semilogy(t_bdf1,err_bdf1(2,:),'b'); hold on; 
%         semilogy(t_rk1,err_rk1(2,:),'b--');
%         semilogy(t_glm1,err_glm1(2,:),'b:');
%         xlabel('t');
%         ylabel('error in v_C');
%         axis tight;
%         subplot(2,3,6)
%         semilogy(t_bdf1,err_bdf1(3,:),'b'); hold on; 
%         semilogy(t_rk1,err_rk1(3,:),'b--');
%         semilogy(t_glm1,err_glm1(3,:),'b:');
%         xlabel('t');
%         ylabel('error in v_R');
%         axis tight; 
        
