% Copyright (c) Tobias Bruell.  All rights reserved.
function solve_dae(Efun,Afun,ffun,tb,t0,tf,xhat,h,tol,disc_method,solfun)

global xxx_Efun;
global xxx_Afun;
global xxx_ffun;
global xxx_h;
global xxx_discretization_method;
xxx_Efun = Efun; xxx_Afun = Afun; xxx_ffun = ffun;
xxx_h = h; xxx_discretization_method = disc_method;

kb = floor(tb/h);
kf = ceil(tf/h);
k0 = round(t0/h);

tic
[x,isunique,r_f,h_f,mu_f,r_b,h_b,mu_b]=...
 solve_dds('solve_dae_Efun','solve_dae_Afun','solve_dae_ffun',kb,k0,kf,xhat,tol);
 toc

r_h_mu___f = [r_f , nan , h_f , nan , mu_f ];
r_h_mu___b = [r_b , nan , h_b , nan , mu_b ];

 isunique

if( nargin > 10 )
   realx = [];
   nsum = 0;
   maxerr = 0;
   for k=kb:kf
      realx = [realx,feval(solfun,k*h)];
      nsum = nsum + norm( x(:,k-kb+1) - realx(:,k-kb+1) );
      if( maxerr < norm( x(:,k-kb+1) - realx(:,k-kb+1) ) )
         maxerr = norm( x(:,k-kb+1) - realx(:,k-kb+1) );
      end
   end
   avg_sol_diff = nsum / (kf-kb)
   max_sol_diff = maxerr

   realx = [];
   fineh = ((tf-tb)/10000);
   for finet=tb:fineh:tf
      realx = [realx,feval(solfun,finet)];
   end

   hold off;
   plot( (kb:kf)*h, x(1,:),     '--', 'Color', [0 0 0]);
   hold on;
   plot( (tb:fineh:tf), realx, '-', 'Color', [0 0 0]);
   plot( (kb:kf)*h, x(2:size(x,1),:),     '--', 'Color', [0 0 0]);
   hold off;
   leg = {};
%   for i = 1:size(x,1)
%      leg{i} = ['computed x',num2str(i)];
%      leg{i+size(x,1)} = ['  real   x',num2str(i)];
%   end
      leg{1} = 'computed solution';
      leg{2} = 'real solution';
      leg
   legend(leg);
   xlabel(['h = ',num2str(h)]);
else
   plot( (kb:kf)*h, x )
end

nsum = 0;
for k=kb:kf-1
   nsum = nsum + norm( feval('solve_dae_Efun',k) * x(:,k+1-kb+1) ...
                     - feval('solve_dae_Afun',k) * x(:,k  -kb+1) ...
                     - feval('solve_dae_ffun',k) );
end
avg_local_error = nsum / (kf-kb)
