% Runme2

 path(path,'~/work/codes/program/Examples');

choice = 1;choicemax = 10;
while( choice < choicemax+1 )
    choice=menu('Beispielaufgabe ',...
        'Example 1 - original EoM',...
        'Truck',...
        'Sand',...
        'Circuit',...
        'Heaviside',...
        'tanh',...
        'tanh_rot',...
        'Example 1 - rot',...
        'Example 2 - linear',...
        'Sand - transformed',...
        'Ende');


    if(choice<choicemax+1)

        % Example 1
        if choice==1
            clear t h y yp ye y0 y1 t0 err ye_2 T y_i err_2 err_i ATOL RTOL i a;
            figure(1);clf; figure(2);clf; figure(3);clf; figure(4);clf;figure(5);clf;

            funct = 'example1';
            funct_exact = 'example1_exact';
            funct_para ='example1_para';
            funct_d1 ='example1_d1';
            funct_d2 = 'example1_d2';
            funct_d3 = 'example1_d3';
            funct_start = 'example1_start';

            a = 0.5; 
            h = 0.001;
            
           RTOL = [ 1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4];

           ATOL = [ 1e-4
                    1e-4
                    1e-2
                    1e-4
                    1e-4
                    1e-1
                    1e-4
                    1e-4
                    1e-3];
            

        end;  %choice


        % Example Truck
        if choice==2
            clear t h y ye y0 y1 t0 err ye_2 T y_i err_2 ATOL RTOL ;
            figure(1);clf

            funct = 'truck';
            funct_d1 = 'truck_d1';
            funct_d2 = 'truck_d2';
            funct_d3 = 'truck_d3';
            funct_exact = 'truck_exact';
            funct_para ='truck_para';

            % Truck :
            a = 20.0;
            
            ATOL=1e-4*ones(n,1);
            ATOL(n)=1e-2;
            ATOL(n+1:n+n)=1e-3*ones(n,1);

            RTOL=1e-3*ones(n,1);
            RTOL(n)=1e-0;
            RTOL(n+1:n+n)=1e-2*ones(n,1);

        end;  %choice


        % Example Sand
        if choice==3
            clear t h y yp ye y0 y1 t0 err ye_2 T y_i err_2 err_i ATOL RTOL i a;
            figure(1);clf;figure(2);clf; figure(3);clf; figure(4);clf;figure(5);clf;

            funct = 'sand';
            funct_d1 = 'sand_d1';
            funct_d2 = 'sand_d2';
            funct_d3 = 'sand_d3';
            funct_exact = 'sand_exact';
            funct_para ='sand_para';
            funct_start = 'sand_start';

            a = 0.5;
            h = 0.001;
            
%             RTOL = [1e-4
%                     1e-4
%                     1e-4
%                     1e-4
%                     1e-4
%                     1e-4
%                     1e-4
%                     1e-4
%                     1e-4];
% 
%             ATOL = [1e-4
%                     1e-4
%                     1e-2
%                     1e-4
%                     1e-4
%                     1e-1
%                     1e-4
%                     1e-4
%                     1e-3];
            
            
        end;  %choice


        % Circuit
        if choice==4
            clear t h y yp ye y0 y1 t0 err ye_2 T y_i err_2 err_i ATOL RTOL;
            figure(1);clf

            funct = 'circuit';
            funct_d1 = 'circuit_d1';
            funct_d2 = 'circuit_d2';
            funct_d3 = 'circuit_d3';
            funct_exact = 'circuit_exact';
            funct_para ='circuit_para';
            funct_start = 'circuit_start';

            a = 7;
            h = 0.01;
            
        end;  %choice
        
        
        % Heaviside
        if choice==5
            clear t h y yp ye y0 y1 t0 err ye_2 T y_i err_2 err_i ATOL RTOL;
            figure(1);clf

            funct = 'heavyside';
            funct_d1 = 'heavyside_d1';
            funct_d2 = 'heavyside_d2';
            funct_d3 = 'heavyside_d3';
            funct_exact = 'heavyside_exact';
            funct_para ='heavyside_para';
            funct_start = 'heavyside_start';

            a = 2;
            h = 0.001;
            
            RTOL=[1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5];
              
            ATOL=[1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5];      
            
        end;  %choice
        
        % tanh
        if choice==6
            clear t h y yp ye y0 y1 t0 err ye_2 T y_i err_2 err_i ATOL RTOL;
            figure(1);clf

            funct = 'ftanh';
            funct_d1 = 'tanh_d1';
            funct_d2 = 'tanh_d2';
            funct_d3 = 'tanh_d3';
            funct_exact = 'tanh_exact';
            funct_para ='tanh_para';
            funct_start = 'tanh_start';

            a = 2;
            h = 0.001;
            
            RTOL=[1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5];
              
            ATOL=[1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5];      
            
        end;  %choice
        
        
        % tanh_rot
        if choice==7
            clear t h y yp ye y0 y1 t0 err ye_2 T y_i err_2 err_i ATOL RTOL;
            figure(1);clf

            funct = 'tanh_rot';
            funct_d1 = 'tanhrot_d1';
            funct_d2 = 'tanhrot_d2';
            funct_d3 = 'tanhrot_d3';
            funct_exact = 'tanhrot_exact';
            funct_para ='tanhrot_para';
            funct_start = 'tanhrot_start';

            a = 2;
            h = 0.001;
            
            RTOL=[1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5];
              
            ATOL=[1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5
                  1e-5];      
            
        end;  %choice
        
      % Example 1- rot
        if choice==8
            clear t h y yp ye y0 y1 t0 err ye_2 T y_i err_2 err_i ATOL RTOL i a;
            figure(1);clf; figure(2);clf; figure(3);clf; figure(4);clf;figure(5);clf;

            funct = 'example1_rot';
            funct_exact = 'example1_exact';
            funct_para ='example1_para';
            funct_d1 ='example1_rot_d1';
            funct_d2 = 'example1_rot_d2';
            funct_d3 = 'example1_rot_d3';
            funct_start = 'example1_start';

            a = 0.3; 
            h = 0.001;
            
           RTOL = [ 1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4];

           ATOL = [ 1e-4
                    1e-4
                    1e-2
                    1e-4
                    1e-4
                    1e-1
                    1e-4
                    1e-4
                    1e-3];
            

        end;  %choice        
        
        % Example2 linear
        if choice==9
            clear t h y yp ye y0 y1 t0 err ye_2 T y_i err_2 err_i ATOL RTOL;
            figure(1);clf;figure(2);clf;figure(3);clf;

            funct = 'ex2Lin';
            funct_d1 = 'ex2Lin_d1';
            funct_d2 = 'ex2Lin_d2';
            funct_d3 = 'ex2Lin_d3';
            funct_exact = 'ex2Lin_exact';
            funct_para ='ex2Lin_para';
            %funct_start = 'ex2Lin_start';

            a = 2;
            h = 0.001;
            
            RTOL=[1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4];
              
            ATOL=[1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4];      
            
        end;  %choice     
        
        
       % Example Sand-transformed
        if choice==10
            clear t h y yp ye y0 y1 t0 err ye_2 T y_i err_2 err_i ATOL RTOL i a;
            figure(1);clf;figure(2);clf; figure(3);clf; figure(4);clf;figure(5);clf;

            funct = 'sand_transf';
            funct_d1 = 'sand_transf_d1';
            funct_d2 = 'sand_transf_d2';
            funct_d3 = 'sand_transf_d3';
            funct_exact = 'sand_transf_exact';
            funct_para ='sand_transf_para';
            %funct_start = 'sand_start';

            a = 0.5;
            h = 0.001;
            
            RTOL = [1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4
                    1e-4];

            ATOL = [1e-4
                    1e-4
                    1e-2
                    1e-4
                    1e-4
                    1e-1
                    1e-4
                    1e-4
                    1e-3];
            
            
        end;  %choice
        
        
        
        

        method_selection = menu('method',...
            'BDF method',...
            'RK method',...
            'GLM method');


        stepsize_selection = menu('stepsize',...
            'variable stepsizes',...
            'constant stepsizes');
        step = 1;
        if stepsize_selection==1
            step = 1;
        end
        if stepsize_selection==2
            step = 2;
        end


        if choice == 3  % Sand
            if method_selection ==1 % bdf
                RTOL=[1e-4
                      1e-4
                      1e-4
                      1e-4
                      1e-4
                      1e-4];

                ATOL=[ 1e-4
                       1e-4
                       1e-3
                       1e-3
                       1e-3
                       1e-3];
            end
            if method_selection ==2 % RK
               RTOL=[ 1e-4
                      1e-4
                      1e-4
                      1e-4
                      1e-4
                      1e-0];

                ATOL=[ 1e-4
                       1e-4
                       1e-4
                       1e-4
                       1e-4
                       1e-4];
            end
            if method_selection ==3 % GLM
                ATOL = 1e-3*ones(12,1);
                RTOL = ATOL;
            end
        end

        if choice == 4  % Circuit
            if method_selection ==1 % bdf
             RTOL=[1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4];
              
            ATOL=[1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4];     
                  
                           
            else
            RTOL=[1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4];
               
            ATOL=[1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4]; 
            end
            
        end
        

        if method_selection == 1
            [t0,y0,y1]=feval(funct_para,'bdf');
            n=length(y0);
            nL = n-length(y1);
            %time = clock;
            [t,h,y] = bdfk2_dae(funct,funct_d1,funct_d2,funct_d3,nL,y0,y1,t0,a,10,2,RTOL,ATOL,h,step);
            %HOBDF(F,dF,ddF,N,T,TOUT,X,XPRIME,CVAL,IPAR,RPAR,IWORK,LIW,RWORK,LRW,RTOL,ATOL,METHOD,INFO,IWARN,IERR);
            %etime(clock, time)
            t_bdf2=t;
            h_bdf2=h;
            y_bdf2=y;
            clear err_bdf2;
            
            for i=1:length(t_bdf2)
                ye(:,i)=feval(funct_exact,t_bdf2(i));
            end
            for i=1:3
                err_bdf2(i,:)=abs(y_bdf2(i,:)-ye(i,:));
            end
            
            err = err_bdf2;
        end
        if method_selection == 2
            [t0,y0,y1]=feval(funct_para,'rk');
            %n=length(y0);
            %time = clock;
            [t,h,y,y_i,T] = RK2dae(funct,funct_d1,funct_d2,funct_d3,y0,y1,t0,a,10,2,RTOL,ATOL,h,'Radau',step);
            %etime(clock, time)
            t_rk2=t;
            h_rk2=h;
            y_rk2=y;
            clear err_rk2;
                        
            for i=1:length(t_rk2)
                ye(:,i)=feval(funct_exact,t_rk2(i));
            end
            for i=1:3
                err_rk2(i,:)=abs(y_rk2(i,:)-ye(i,:));
            end
            
            err = err_rk2;
            
        end
        if method_selection == 3
            [t0,y0,y1]=feval(funct_para,'glm');
            %n=length(y0);
            p=2;
            %time = clock;
            [t,h,y,yp,y_est]=glm2_dae(funct,funct_start,funct_d1,funct_d2,funct_d3,y0,y1,t0,a,10,p,p,p+1,p+1,RTOL,ATOL,h,step);
            %etime(clock, time)
            t_glm2=t;
            h_glm2=h;
            y_glm2=y;
            clear err_glm2;
            
            for i=1:length(t_glm2)
                ye(:,i)=feval(funct_exact,t_glm2(i));
            end
            for i=1:3
                err_glm2(i,:)=abs(y_glm2(i,:)-ye(i,:));
            end
            
            err = err_glm2;
        end

%         for i=1:length(t)
%             ye(:,i)=feval(funct_exact,t(i));
%         end
%         [m,n]=size(ye);
%         for i=1:m
%             err(i,:)=abs(y(i,:)-ye(i,:));
%         end


        plot_results(t,y,ye,err);
        % output

%         figure(1);
%         subplot(2,2,1)
%         plot(t,ye(1,:),'r'); hold on;
%         plot(t,y(1,:),'b');
%         xlabel('t');
%         ylabel('x');
%         axis tight;
%         subplot(2,2,2)
%         plot(t,ye(2,:),'r'); hold on;
%         plot(t,y(2,:),'b');
%         xlabel('t');
%         ylabel('y');
%         axis tight;
%         subplot(2,3,3)
%         plot(t,ye(3,:),'r'); hold on;
%         plot(t,y(3,:),'b');
%         axis tight;
%         xlabel('t');
%         ylabel('L');


%         subplot(2,2,3)
%         semilogy(t,err(1,:),'b--');
%         xlabel('t');
%         ylabel('x');
%         axis tight;
%         subplot(2,2,4)
%         semilogy(t,err(2,:),'b--');
%         xlabel('t');
%         ylabel('y');
%         axis tight;
%         subplot(2,3,6)
%         semilogy(t,err(3,:),'b--');
%         xlabel('t');
%         ylabel('L');
%         axis tight;

    end

end;  %while


if(1==menu('Grafiken','Grafiken schliessen','Grafiken lassen')), close all; end;

! more bdfk2_dae.m bdf2_step.m RK2dae.m  RK_step.m  glm2_dae.m  glm_external.m glm_internal.m starting.m startProc.m >src.txt

%bdfk1_dae.m bdf1.m   RK1dae.m  RK_step1.m   glm1dae.m  glm1_external.m  glm1_internal.m starting1.m

! mpage -4A src.txt > src_4p.ps


disp(' ');






