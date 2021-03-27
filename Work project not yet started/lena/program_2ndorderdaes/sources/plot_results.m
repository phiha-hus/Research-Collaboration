function plot_results(t,y,ye,err)

[m,n]=size(ye);

for i=1:3

    figure(i);
    subplot(2,1,1)
    plot(t,ye(i,:),'r'); hold on;
    plot(t,y(i,:),'b');
    xlabel('t');
    ylabel('x');
    axis tight;
    
    subplot(2,1,2)
    semilogy(t,err(i,:),'b--');
    xlabel('t');
    ylabel('x');
    axis tight;
end

