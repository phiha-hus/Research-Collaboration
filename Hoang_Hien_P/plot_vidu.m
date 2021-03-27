%% Plot 2 vi du trong bai bao triggering

%% 
figure(1); clf;
semilogy(T_plot,X_plot,'-.');
grid on
title('Plot the solution w.r.t piecewise constant input u')    
xlabel('t (sec.)')
ylabel('x(t)')
legend('x_1','x_2','x_3')

print(gcf, '-dpdf', 'Fig_1a.pdf');


%%

N = length(T_trig);
T_trig2 = T_trig(2:end) - T_trig(1:end-1) ; 

figure(2); clf;
plot(T_trig(1:end-1),T_trig2,'r *'); hold on
grid on
%legend('t\_trig')
xlabel('t(sec.)')
ylabel('T_k = t_{k+1}-t_k')
title('Find all triggering point on [0,10]');
hold off;

print(gcf, '-dpdf', 'Fig_1b.pdf');

%% 
figure(3); clf;

subplot(2,1,1)

for i = 1:N-1
    A = [T_trig(i) T_trig(i)];
    B = [0 T_trig2(i)] ;
    line(A,B,'Color','red','LineStyle','-'); hold on  
    
    % plot a small red circle around the point    
    plot(T_trig(i),T_trig2(i), 'Marker', 'o','MarkerEdgeColor','r','MarkerSize',4); hold on
end

grid on
%legend('t\_trig')
xlabel('t(sec.)')
ylabel('T_k = t_{k+1}-t_k')
title('Find all triggering point on [0,10]');
hold off;

print(gcf, '-dpdf', 'Fig_1c.pdf');

%%
figure(4); clf;

subplot(2,1,1)
semilogy(T_plot,X_plot,'-.');
grid on
title('Plot the solution w.r.t piecewise constant input u')    
xlabel('t (sec.)')
ylabel('x(t)')
legend('x_1','x_2','x_3')
% resizeLegend();


subplot(2,1,2)

for i = 1:N-1
    A = [T_trig(i) T_trig(i)];
    B = [0 T_trig2(i)] ;
    line(A,B,'Color','red','LineStyle','-'); hold on  
    
    % plot a small red circle around the point    
    plot(T_trig(i),T_trig2(i), 'Marker', 'o','MarkerEdgeColor','r','MarkerSize',3); hold on
end

grid on
%legend('t\_trig')
xlabel('t(sec.)')
ylabel('T_k = t_{k+1}-t_k')
title('Find all triggering point on [0,10]');
hold off;

print(gcf, '-dpdf', 'Fig_1d.pdf');

