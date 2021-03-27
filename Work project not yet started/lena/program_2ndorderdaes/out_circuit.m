      figure(1);
        subplot(2,3,1)
        plot(t,ye(1,:),'r'); hold on; 
        plot(t_bdf,y_bdf(1,:),'b'); 
        plot(tk,y_rk(1,:),'b--'); 
        plot(t_glm,y_glm(1,:),'b:'); 
   
        xlabel('t'); 
        ylabel('v_L');
        axis tight;
        subplot(2,3,2)
        plot(t,ye(2,:),'r'); hold on; 
        plot(t_bdf,y_bdf(2,:),'b'); 
        plot(tk,y_rk(2,:),'b--'); 
        plot(t_glm,y_glm(2,:),'b:'); 
        xlabel('t');
        ylabel('v_C');
        axis tight;  
        subplot(2,3,3)
        plot(t,ye(3,:),'r'); hold on; 
        plot(t_bdf,y_bdf(3,:),'b'); 
        plot(tk,y_rk(3,:),'b--'); 
        plot(t_glm,y_glm(3,:),'b:'); 
        xlabel('t');
        ylabel('v_R');
        axis tight;  
  
        
        subplot(2,3,4)
        semilogy(t_bdf,err_bdf(1,:),'b'); hold on; 
        semilogy(tk,err_rk(1,:),'b--');
        semilogy(t_glm,err_glm(1,:),'b:');
        xlabel('t');
        ylabel('error in v_L');
        axis tight;
        subplot(2,3,5)
        semilogy(t_bdf,err_bdf(2,:),'b'); hold on; 
        semilogy(tk,err_rk(2,:),'b--');
        semilogy(t_glm,err_glm(2,:),'b:');
        xlabel('t');
        ylabel('error in v_C');
        axis tight;
        subplot(2,3,6)
        semilogy(t_bdf,err_bdf(3,:),'b'); hold on; 
        semilogy(tk,err_rk(3,:),'b--');
        semilogy(t_glm,err_glm(3,:),'b:');
        xlabel('t');
        ylabel('error in v_R');
        axis tight;
        

        
      
        