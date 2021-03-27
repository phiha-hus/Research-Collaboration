      figure(1);
        title('constant');
        subplot(2,3,1)
        %plot(t,ye(1,:),'r'); hold on; 
        plot(t_glm_c,y_glm_c(1,:),'b--');hold on; 
        plot(t_glm1_c,y_glm1_c(1,:),'b');
        xlabel('t'); 
        ylabel('x');
        axis tight;
        subplot(2,3,2)
        %plot(t,ye(2,:),'r'); hold on;    
        plot(t_glm_c,y_glm_c(2,:),'b--');hold on; 
        plot(t_glm1_c,y_glm1_c(2,:),'b');
        xlabel('t');
        ylabel('y');
        axis tight;  
        subplot(2,3,3)
        %plot(t,ye(3,:),'r'); hold on;    
        plot(t_glm_c,y_glm_c(3,:),'b--');hold on; 
        plot(t_glm1_c,y_glm1_c(3,:),'b');
        xlabel('t');
        ylabel('\lambda');
        axis tight;  
  
        
        subplot(2,3,4)
        semilogy(t_glm_c,err_glm_c(1,:),'b--'); hold on; 
        semilogy(t_glm1_c,err_glm1_c(1,:),'b');
        xlabel('t');
        ylabel('error in x');
        axis tight;
        subplot(2,3,5)
        semilogy(t_glm_c,err_glm_c(2,:),'b--'); hold on; 
        semilogy(t_glm1_c,err_glm1_c(2,:),'b');
        xlabel('t');
        ylabel('error in y');
        axis tight;
        subplot(2,3,6)
        semilogy(t_glm_c,err_glm_c(3,:),'b--'); hold on; 
        semilogy(t_glm1_c,err_glm1_c(3,:),'b');
        xlabel('t');
        ylabel('error in \lambda');
        axis tight;
        

        
        figure(2);
        title('variable');
        
        subplot(2,3,1)
        %plot(t,ye(1,:),'r'); hold on; 
        plot(t_glm_var,y_glm_var(1,:),'b--');hold on; 
        plot(t_glm1_var,y_glm1_var(1,:),'b');
        xlabel('t'); 
        ylabel('x');
        axis tight;
        subplot(2,3,2)
        %plot(t,ye(2,:),'r'); hold on;    
        plot(t_glm_var,y_glm_var(2,:),'b--');hold on; 
        plot(t_glm1_var,y_glm1_var(2,:),'b');
        xlabel('t');
        ylabel('y');
        axis tight;  
        subplot(2,3,3)
        %plot(t,ye(3,:),'r'); hold on;    
        plot(t_glm_var,y_glm_var(3,:),'b--');hold on; 
        plot(t_glm1_var,y_glm1_var(3,:),'b');
        xlabel('t');
        ylabel('\lambda');
        axis tight;  
  
        
        subplot(2,3,4)
        semilogy(t_glm_var,err_glm_var(1,:),'b--'); hold on; 
        semilogy(t_glm1_var,err_glm1_var(1,:),'b');
        xlabel('t');
        ylabel('error in x');
        axis tight;
        subplot(2,3,5)
        semilogy(t_glm_var,err_glm_var(2,:),'b--'); hold on; 
        semilogy(t_glm1_var,err_glm1_var(2,:),'b');
        xlabel('t');
        ylabel('error in y');
        axis tight;
        subplot(2,3,6)
        semilogy(t_glm_var,err_glm_var(3,:),'b--'); hold on; 
        semilogy(t_glm1_var,err_glm1_var(3,:),'b');
        xlabel('t');
        ylabel('error in \lambda');
        axis tight;
        