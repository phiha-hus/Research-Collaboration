

[m,n]=size(ye);

figure(1);
%   subplot(2,3,1)
%     %plot(t,ye(1,:),'r'); hold on;
%     plot(t_bdf2,y_bdf2(1,:),'b');hold on;
%     plot(t_bdf1,y_bdf1(1,:),'b--');
%     plot(t_rk2,y_rk2(1,:),'r');
%     plot(t_rk1,y_rk1(1,:),'r--');
%     plot(t_glm2,y_glm2(1,:),'c');
%     plot(t_glm1,y_glm1(1,:),'c--');   
%     xlabel('t');
%     ylabel('x');
%     axis tight;
%     
%   subplot(2,3,2)
%    % plot(t,ye(2,:),'r'); hold on;
%     plot(t_bdf2,y_bdf2(2,:),'b');hold on;
%     plot(t_bdf1,y_bdf1(2,:),'b--');
%     plot(t_rk2,y_rk2(2,:),'r');
%     plot(t_rk1,y_rk1(2,:),'r--');
%     plot(t_glm2,y_glm2(2,:),'c');
%     plot(t_glm1,y_glm1(2,:),'c--');   
%     xlabel('t');
%     ylabel('x');
%     axis tight;    
%     
%    subplot(2,3,3)
%     %plot(t,ye(3,:),'r'); hold on;
%     plot(t_bdf2,y_bdf2(3,:),'b');hold on;
%     plot(t_bdf1,y_bdf1(3,:),'b--');
%     plot(t_rk2,y_rk2(3,:),'r');
%     plot(t_rk1,y_rk1(3,:),'r--');
%     plot(t_glm2,y_glm2(3,:),'c');
%     plot(t_glm1,y_glm1(3,:),'c--');   
%     xlabel('t');
%     ylabel('x');
%     axis tight;   
    
    subplot(1,3,1)
    semilogy(t_bdf2,err_bdf2(1,:),'b');hold on;
    semilogy(t_bdf1,err_bdf1(1,:),'b--');
    semilogy(t_rk2,err_rk2(1,:),'r');
    semilogy(t_rk1,err_rk1(1,:),'r--');
    %semilogy(t_glm2,err_glm2(1,:),'c');
    %semilogy(t_glm1,err_glm1(1,:),'c--');    
    xlabel('t');
    ylabel('x');
    axis tight;

    subplot(1,3,2)
    semilogy(t_bdf2,err_bdf2(2,:),'b');hold on;
    semilogy(t_bdf1,err_bdf1(2,:),'b--');
    semilogy(t_rk2,err_rk2(2,:),'r');
    semilogy(t_rk1,err_rk1(2,:),'r--');
    %semilogy(t_glm2,err_glm2(2,:),'c');
    %semilogy(t_glm1,err_glm1(2,:),'c--');    
    xlabel('t');
    ylabel('x');
    axis tight;
    
    subplot(1,3,3)
    semilogy(t_bdf2,err_bdf2(3,:),'b');hold on;
    semilogy(t_bdf1,err_bdf1(3,:),'b--');
    semilogy(t_rk2,err_rk2(3,:),'r');
    semilogy(t_rk1,err_rk1(3,:),'r--');
    %semilogy(t_glm2,err_glm2(3,:),'c');
    %semilogy(t_glm1,err_glm1(3,:),'c--');    
    xlabel('t');
    ylabel('x');
    axis tight;
