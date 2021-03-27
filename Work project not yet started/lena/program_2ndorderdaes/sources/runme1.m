% runme 1
 path(path,'~/work/codes/program/Examples');
 %path(path,'~/program/1st order');
 %path(path,'/1st order');


choice=1;choicemax=9;
while(choice<choicemax+1)
    choice=menu('Beispielaufgabe ',...
        'Example 1 - original EoM',...
        'Truck',...
        'Sand',...
        'Sand - strangenessfree',...
        'Circuit',...
        'Steffen',...
        'Heaviside',...
        'Heaviside - strangenessfree',...
        'Example 1 - rotate',...
        'Ende');
    
    if(choice<choicemax+1)
          
        % Example 1 
        if choice==1 
            clear t h y yp ye y0 y1 t0 err ye_2 T y_i err_2 err_i;
            figure(1);clf;
            
            funct = 'example1_1';
            funct_d1 ='example1_1_d1';
            funct_d2 = 'example1_1_d2';
            funct_exact = 'example1_1_exact';
            funct_para ='example1_1_para';
            funct_start = 'example1_1_start';
            
            a = 0.5;
            N = 10; 
            h = 0.001;
            
            RTOL=[1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-0
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4];
            
            ATOL=[1e-4
                  1e-4
                  1e-2 % 1
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-3
                  1e-1
                  1e-4
                  1e-2
                  1e-4
                  1e-4
                  1e-4];   
            
           
        end;  %choice
       
        
        % Example Truck 
        if choice==2
            clear t h y ye y0 y1 t0 err ye_2 T y_i err_2;
            figure(1);clf
  
            funct = 'example_truck_1';
            funct_exact = 'example_truck_exact_1';
            funct_para ='example_truck_para_1';
       
            % Truck :
            a = 20.0; 
  
        end;  %choice
        
        
        % Example Sand 
        if choice==3
            clear t h y ye y0 y1 t0 err ye_2 T y_i err_2 err_i;
            figure(1);clf;
            
            funct = 'sand_1';
            funct_d1 ='sand_1_d1';
            funct_d2 = 'sand_1_d2';
            funct_exact = 'sand_1_exact';
            funct_para ='sand_1_para';
            funct_start = 'sand_1_start'; 
            
            
            a = 0.5;
            N = 10;
            h = 0.001;
            
            RTOL=[1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-0
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4];
            
            ATOL=[1e-4
                  1e-4
                  1e-2 % 1
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-3
                  1e-1
                  1e-4
                  1e-2
                  1e-4
                  1e-4
                  1e-4]; 
        end;  %choice
        
        
        % Example Sand - strangnessfree
        if choice==4
            clear t h y ye y0 y1 t0 err ye_2 T y_i err_2;
            figure(1);clf
            
            funct = 'example_sand_1_strf';
            funct_exact = 'example_sand_exact_1_strf';
            funct_para ='example_sand_para_1_strf';
            
            a = 0.5; 
            N = 10;

        end;  %choice
        
        
        % Circuit
        if choice==5
            clear t h y ye y0 y1 t0 err L1  dim  y_2 L2 ;
            figure(1);clf;
           
            funct = 'circuit_1';
            funct_exact = 'circuit_exact_1';
            funct_para ='circuit_para_1';
            funct_d1 ='circuit_1_d1';
            funct_d2 = 'circuit_1_d2';
            funct_start = 'circuit_1_start';
            
         
            a = 7;
            N = 10;
            h = 0.01;
            
            RTOL=[1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
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
                  1e-3
                  1e-3
                  1e-2
                  1e-4
                  1e-4
                  1e-3
                  1e-3
                  1e-2
                  1e-2
                  1e-2
                  1e-3
                  1e-2
                  1e-2
                  1e-2];
      
        end % choice
        
        % Steffen
        if choice==6
            clear t h y ye y0 y1 t0 err L1  dim  y_2 L2 ;
            figure(1);clf;
           
            funct = 'Steffen_Exp';
            funct_exact = 'Steffen_Exp_exact';
            funct_para ='Steffen_Exp_para';
            funct_d1 ='Steffen_Exp_d1';
            funct_d2 = 'Steffen_Exp_d2';
            funct_start = 'Steffen_Exp_start';
            
         
            a = 0.5;
            N = 10;
            h = 0.015;
            
            RTOL=[1e-4
                  1e-4
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
                  1e-3
                  1e-3
                  1e-2
                  1e-4
                  1e-4
                  1e-3
                  1e-3
                  1e-2];
      
        end % choice     
        
        
        % Heaviside
        if choice==7
            clear t h y ye y0 y1 t0 err L1  dim  y_2 L2 ;
            figure(1);clf;
           
            funct = 'heaviside_1';
            funct_exact = 'heaviside_exact_1';
            funct_para ='heaviside_para_1';
            funct_d1 ='heaviside_1_d1';
            funct_d2 = 'heaviside_1_d2';
            funct_start = 'heaviside_1_start';
            
         
            a = 2;
            N = 10;
            h = 0.01;
            
            RTOL=[1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
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
                  1e-3
                  1e-3
                  1e-2
                  1e-4
                  1e-4
                  1e-3
                  1e-3
                  1e-2
                  1e-2
                  1e-2
                  1e-3
                  1e-2
                  1e-2
                  1e-2];
      
        end % choice        % Circuit
        if choice==8
            clear t h y ye y0 y1 t0 err L1  dim  y_2 L2 ;
            figure(1);clf;
           
            funct = 'heaviside_sf_1';
            funct_exact = 'heaviside_sf_exact_1';
            funct_para ='heaviside_sf_para_1';
            funct_d1 ='heaviside_sf_1_d1';
            funct_d2 = 'heaviside_sf_1_d2';
            funct_start = 'heaviside_sf_1_start';
            
         
            a = 2;
            N = 10;
            h = 0.01;
            
            RTOL=[1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
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
                  1e-3
                  1e-3
                  1e-2
                  1e-4
                  1e-4
                  1e-3
                  1e-3
                  1e-2
                  1e-2
                  1e-2
                  1e-3
                  1e-2
                  1e-2
                  1e-2];
      
        end % choice

        % Example 1 - rotate 
        if choice==9 
            clear t h y yp ye y0 y1 t0 err ye_2 T y_i err_2 err_i;
            figure(1);clf; figure(2);clf; figure(3);clf;
            
            funct = 'example1_rot_1';
            funct_d1 ='example1_rot_1_d1';
            funct_d2 = 'example1_rot_1_d2';
            funct_exact = 'example1_1_exact';
            funct_para ='example1_1_para';
           % funct_start = 'example1_1_start';
            
            a = 0.3;
            N = 10; 
            h = 0.001;
            
            RTOL=[1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-0
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4];
            
            ATOL=[1e-4
                  1e-4
                  1e-2 % 1
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-4
                  1e-3
                  1e-1
                  1e-4
                  1e-2
                  1e-4
                  1e-4
                  1e-4];   
            
           
        end;  %choice        
       method_selection = menu('method',...
            'BDF method',...
            'RK method',...
            'GLM method');
        
       
        stepsize_selection = menu('stepsize',...
             'variable stepsizes',...
             'constant stepsizes');
        if stepsize_selection==1
            step = 1;
        end
        if stepsize_selection==2
            step = 2;
        end
        
        
        if method_selection ==1
            [t0,y0]=feval(funct_para);
            n=length(y0);
            [t,h,y,y_est,y_i]=bdfk1_dae(funct,funct_d1,funct_d2,y0,t0,a,N,RTOL,ATOL,h,step,3); 

            y_bdf1=y;
            t_bdf1=t;          
            clear err_bdf1;
            
            for i=1:length(t_bdf1)
                ye(:,i)=feval(funct_exact,t_bdf1(i));
            end
            for i=1:3
                err_bdf1(i,:)=abs(y_bdf1(i,:)-ye(i,:));
            end
            
            err=err_bdf1;
            
        end
        if method_selection ==2
            [t0,y0,y1]=feval(funct_para);
            n=length(y0);
            [t,h,y,y_i,T]=RK1dae(funct,funct_d1,funct_d2,y0,t0,a,N,2,RTOL,ATOL,h,'Radau',step);

            y_rk1=y;
            t_rk1=t;
            clear err_rk1;
            
            for i=1:length(t_rk1)
                ye(:,i)=feval(funct_exact,t_rk1(i));
            end
            for i=1:3
                err_rk1(i,:)=abs(y_rk1(i,:)-ye(i,:));
            end
            
            err = err_rk1;
        end
        if method_selection ==3
            [t0,y0]=feval(funct_para);
            n=length(y0);
            [t,h,y,yp]=glm1dae(funct,funct_start,funct_d1,funct_d2,y0,t0,a,N,2,2,3,3,RTOL,ATOL,h,step);

            y_glm1=y;
            t_glm1=t;
            clear err_glm1;
            
            for i=1:length(t_glm1)
                ye(:,i)=feval(funct_exact,t_glm1(i));
            end
            for i=1:3
                err_glm1(i,:)=abs(y_glm1(i,:)-ye(i,:));
            end
            
            err=err_glm1;
        end

        
        figure(1);
        subplot(2,3,1)
        plot(t,ye(1,:),'r'); hold on; 
        plot(t,y(1,:),'b*');
        xlabel('t'); 
        ylabel('x');
        axis tight;
        subplot(2,3,2)
        plot(t,ye(2,:),'r'); hold on;    
        plot(t,y(2,:),'b*');
        xlabel('t');
        ylabel('y');
        axis tight;  
        subplot(2,3,3)   
        plot(t,ye(3,:),'r'); hold on; 
        plot(t,y(3,:),'b*');
        axis tight;
        xlabel('t');
        ylabel('L');
        
        
        subplot(2,3,4)
        semilogy(t,err(1,:),'b--');
        xlabel('t');
        ylabel('x');
        axis tight;
        subplot(2,3,5)
        semilogy(t,err(2,:),'b--');
        xlabel('t');
        ylabel('y');
        axis tight;
        subplot(2,3,6)
        semilogy(t,err(3,:),'b--');
        xlabel('t');
        ylabel('L');
        axis tight;
        
        
    end
    
end;  %while


if(1==menu('Grafiken','Grafiken schliessen','Grafiken lassen')), close all; end;

! more bdf1.m bdfk2_dae_1.m bdfk2_dae.m bdfk_coef_1.m bdfk_coef_a.m bdfk_coef.m divDiff.m example1.m EXAMPLE1.M example1_RK.m example_sand_exact.m example_Sand.m example_sand_para.m extrapol.m faculty.m gaussZerl.m interpol.m myNorm.m Newton.m Newton_RK.m RK2dae.m RK_coef.m RK_step.m runme_bdf.m runme_RK.m >src.txt   

! mpage -4A src.txt > src_4p.ps
! mpage -4A res.txt > res_4p.ps


disp(' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%