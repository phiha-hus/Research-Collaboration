% Copyright (c) Tobias Bruell.  All rights reserved.
%
% solve_dds Compute a solution of a discrete-time descriptor system.
%     solve_dds(Efun,Afun,ffun,kb,k0,kf,xhat,tol) tries to compute the
%     iterates x_kb , ... , x_kf of a solution of one of the discrete-time
%     descriptor systems
%
%      (1)  E_k x_{k+1} = A_k x_{k} + f_k    for      kb <= k <= infty
%      (2)  E_k x_{k+1} = A_k x_{k} + f_k    for  -infty <= k <= kf
%      (3)  E_k x_{k+1} = A_k x_{k} + f_k    for  -infty <= k <= infty
%
%     along with         xhat = x_{k0},
%     where the matrices E_k and A_k and the vector f_k are obtained by
%     evaluation the function Efun, Afun and ffun, respectively, i.e.
%         E_k = feval( Efun, k ), A_k = feval( Afun, k ),
%         f_k = feval( ffun, k ).
%     If no such solution exists, the right hand side f_k is changed so
%     that a solution exists. If there are multiple solutions, one
%     solution is selected. If the initial value 'xhat' is inconsistent,
%     the consistent initial value is choosen which is closest to
%     'xhat' in the 2-norm.
%     If the constant rank assumptions form [1] are not satisfied an
%     error message is issued and the algorithm is aborted, although
%     there might still exist a solution.
%     If kb == k0 and k0  < kf equation (1) is considered.
%     If kb  < k0 and k0 == kf equation (2) is considered.
%     If kb  < k0 and k0  < kf equation (3) is considered.
%     If kb == k0 and k0 == kf no equation  is considered and no real
%     computations are performed.
%
%     Input parameters:
%
%        Efun  :  A function handle or name of a function that takes
%                 one input parameter k and returns the matrix E_k.
%
%        Afun  :  A function handle or name of a function that takes
%                 one input parameter k and returns the matrix A_k.
%
%        ffun  :  A function handle or name of a function that takes
%                 one input parameter k and returns the vector f_k.
%
%        kb    :  The index of the first iterate to compute.
%        k0    :  The index, where to apply the initial condition.
%                 k0 has to satisfy kb <= k0 <= kf.
%        kf    :  The index of the last iterate to compute.
%
%        xhat  :  The desired value of the solution at iterate x_k0.
%
%        tol   :  Tolerance, below which a singular value is considered
%                 zero.
%
%     Output parameters
%
%        x     :  The solution iterates x_kb , ... , x_kf side by side
%                 in a matrix. Thus, x(:,k-kb+1) represents the x_k
%                 iterate of the solution. x(:,1) represents the x_kb
%                 iterate and x(kf-kb+1) represents the x_kf iterate.
%
%        isunique : This is set to 1 if there is only one unique
%                 solution and to 0 otherwise.
%
%        r_f   :  The sequence of the r_E as in [1]. Note that r_f(i)
%                 corresponds to r_{E,i-1} in the forward-only
%                 reduction process as described in [1].
%
%        h_f   :  The sequence of the h_f as in [1]. Note that h_f(i)
%                 corresponds to h_{f,i-1} in the forward-only
%                 reduction process as described in [1].
%
%        mu_f  :  The forward strangeness index of the system as
%                 defined in [1].
%
%        r_b   :  The sequence of the r_A as in [1]. Note that r_b(i)
%                 corresponds to r_{A,i-1} in the backward-only
%                 reduction process as described in [1].
%
%        h_b   :  The sequence of the h_b as in [1]. Note that h_b(i)
%                 corresponds to h_{b,i-1} in the backward-only
%                 reduction process as described in [1].
%
%        mu_b  :  The backward strangeness index of the system as
%                 defined in [1].
%
%     Reference:
%        [1]    Bruell, T.
%               Linear discrete-time descriptor systems;
%               Diploma-Thesis (2007);
%               http://www.math.tu-berlin.de/~bruell
%
function [x,isunique,r_f,h_f,mu_f,r_b,h_b,mu_b]=...
            solve_dds(Efun,Afun,ffun,kb,k0,kf,xhat,tol)

r_f=[]; h_f=[]; %(*@\label{cl_det_fb_begin}@*)
mu_f = 0;
r_b=[]; h_b=[];
mu_b = 0;

dimmat = feval( Efun, k0 );
[m,n] = size(dimmat);
clear dimmat;

do_forward  = ( k0 < kf );
do_backward = ( k0 > kb ); %(*@\label{cl_det_fb_end}@*)

if( do_forward ) %(*@\label{cl_det_findex_begin}@*)
   Ek_f = {};
   Ak_f = {};
   Qk_f = {};
   fk_f = {};

   kact = k0;

   % determine the forward index
   haveindex = 0;
   while( ~haveindex )
      [Ek_f,Ak_f,Qk_f,fk_f,r_f,h_f]=advance_inflated_system...
      (Efun,Afun,ffun,m,n,Ek_f,Ak_f,Qk_f,fk_f,kact,mu_f,r_f,h_f,tol,0);

      if( ( mu_f >= 1 && r_f(mu_f+1) == r_f(mu_f) ) )
         haveindex = 1;
      else
         mu_f = mu_f + 1;
      end
   end

   r_mu_f = r_f(mu_f+1);
   h_mu_f = h_f(mu_f+1);
else
   r_mu_f = -1; h_mu_f = -1;
end %(*@\label{cl_det_findex_end}@*)

if( do_backward ) %(*@\label{cl_det_bindex_begin}@*)
   Ek_b = {};
   Ak_b = {};
   Qk_b = {};
   fk_b = {};

   kact = k0-1;

   % determine the backward index
   haveindex = 0;
   while( ~haveindex )
      [Ek_b,Ak_b,Qk_b,fk_b,r_b,h_b]=advance_inflated_system...
      (Afun,Efun,ffun,m,n,Ek_b,Ak_b,Qk_b,fk_b,kact,mu_b,r_b,h_b,tol,1);

      if( ( mu_b >= 1 && r_b(mu_b+1) == r_b(mu_b) ) )
         haveindex = 1;
      else
         mu_b = mu_b + 1;
      end
   end

   r_mu_b = r_b(mu_b+1);
   h_mu_b = h_b(mu_b+1);
else
   r_mu_b = -1; h_mu_b = -1;
end %(*@\label{cl_det_bindex_end}@*)

% find all constraints that the initial condition has to fulfill
constraintA = zeros(0,n);%(*@\label{cl_det_constraints_begin}@*)
constraintf = zeros(0,1);
if( do_forward )%(*@\label{cl_det_constraintsadd_begin}@*)
   constraintA = [constraintA;-Ak_f{1}(r_mu_f+1:r_mu_f+h_mu_f,:)...
                                                        *Qk_f{1}'];
   constraintf = [constraintf; fk_f{1}(r_mu_f+1:r_mu_f+h_mu_f)];
end
if( do_backward )
   constraintA = [constraintA; Ak_b{1}(r_mu_b+1:r_mu_b+h_mu_b,:)...
                                                        *Qk_b{1}'];
   constraintf = [constraintf; fk_b{1}(r_mu_b+1:r_mu_b+h_mu_b)];
end%(*@\label{cl_det_constraintsadd_end}@*)

% find the consistent initial condition that is closest to the given
% one
[U,S,V] = svd(constraintA);
cnum = rank(S); % compute the number of constraints
cnum

ytilde = zeros(n,1);
xtilde = V' * xhat;
ftilde = U' * constraintf;

ytilde(1:cnum)   = S(1:cnum,1:cnum) \ ftilde(1:cnum);
ytilde(cnum+1:n) = xtilde(cnum+1:n);

xhat = V * ytilde;%(*@\label{cl_det_constraints_end}@*)

if( do_forward ) %(*@\label{cl_fsol_begin}@*)
   kact = k0;

   xtt_old = (Qk_f{1})' * xhat;

   % start solver
   for i=1:(kf-k0)
      xtt = zeros(n,1);

      xtt(1:h_mu_f) = -Ak_f{2}(r_mu_f+1:r_mu_f+h_mu_f,1:h_mu_f)\...
                       fk_f{2}( r_mu_f+1 : (r_mu_f+h_mu_f) );

      xtt(h_mu_f+1:h_mu_f+r_mu_f) = ...
            Ek_f{1}(1:r_mu_f,h_mu_f+1:h_mu_f+r_mu_f)\(...
                     Ak_f{1}(1:r_mu_f,:) * xtt_old ...
                   + fk_f{1}(1:r_mu_f) ...
                   - Ek_f{1}(1:r_mu_f,:) * xtt );

      % choose a solution
      xtt(h_mu_f+r_mu_f+1:n) = zeros(n-r_mu_f-h_mu_f,1);

      x(1:n, k0-kb+i+1) = (Qk_f{2}) * xtt;
      xtt_old = xtt;

      % proceed one step
      kact = kact+1;
      for i=1:prod(size(Ek_f))-1
         Ek_f{i}=Ek_f{i+1};
         Ak_f{i}=Ak_f{i+1};
         fk_f{i}=fk_f{i+1};
         Qk_f{i}=Qk_f{i+1};
      end
      [Ek_f,Ak_f,Qk_f,fk_f,r_f,h_f]=advance_inflated_system...
      (Efun,Afun,ffun,m,n,Ek_f,Ak_f,Qk_f,fk_f,kact,mu_f,r_f,h_f,tol,0);
   end

   % compute the _real_ forward strangeness index (i.e. as in [1])
   for i=length(r_f)-1:-1:1
      if( r_f(i) == r_f(i+1) )
         mu_f = i-1;
      end
   end
end %(*@\label{cl_fsol_end}@*)

if( do_backward ) %(*@\label{cl_bsol_begin}@*)
   kact = k0-1;

   xtt_old = (Qk_b{1})' * xhat;

   % start solver
   for i=0:(k0-kb-1)
      xtt = zeros(n,1);%(*@\label{cl_comp_sol_begin}@*)

      xtt(1:h_mu_b) =  Ak_b{2}(r_mu_b+1:r_mu_b+h_mu_b,1:h_mu_b)\...
                       fk_b{2}( r_mu_b+1 : (r_mu_b+h_mu_b) );

      xtt(h_mu_b+1:h_mu_b+r_mu_b) = ...
            Ek_b{1}(1:r_mu_b,h_mu_b+1:h_mu_b+r_mu_b)\( ...  %(*@\label{cl_error1}@*)
                     Ak_b{1}(1:r_mu_b,:) * xtt_old ...
                   - fk_b{1}(1:r_mu_b) ...
                   - Ek_b{1}(1:r_mu_b,:) * xtt );

      % choose a solution
      xtt(h_mu_b+r_mu_b+1:n) = zeros(n-r_mu_b-h_mu_b,1);

      x(1:n, k0-kb-i) =  (Qk_b{2}) * xtt;
      xtt_old = xtt;%(*@\label{cl_comp_sol_end}@*)

      % proceed one step
      kact = kact-1;
      for i=1:prod(size(Ek_b))-1
         Ek_b{i}=Ek_b{i+1};
         Ak_b{i}=Ak_b{i+1};
         fk_b{i}=fk_b{i+1};
         Qk_b{i}=Qk_b{i+1};
      end
      [Ek_b,Ak_b,Qk_b,fk_b,r_b,h_b]=advance_inflated_system...
      (Afun,Efun,ffun,m,n,Ek_b,Ak_b,Qk_b,fk_b,kact,mu_b,r_b,h_b,tol,1);
   end

   % compute the _real_ backward strangeness index (i.e. as in [1])
   for i=length(r_b)-1:-1:1
      if( r_b(i) == r_b(i+1) )
         mu_b = i-1;
      end
   end
end %(*@\label{cl_bsol_end}@*)

if( ( ~do_forward  || r_mu_f + h_mu_f == n ) && ...  %(*@\label{cl_last_begin}@*)
    ( ~do_backward || r_mu_b + h_mu_b == n ) )
   isunique = 1;
else
   isunique = 0;
end

x(:,k0-kb+1) = xhat; %(*@\label{cl_last_end}@*)



function [Ek,Ak,Qk,fk,r,h]=advance_inflated_system...
              (Efun,Afun,ffun,m,n,Ek,Ak,Qk,fk,kact,mu,r,h,tol,backward)

if( backward )
   signed_mu = -mu;
else
   signed_mu =  mu;
end

Ek{mu+1} = feval( Efun, kact+signed_mu );%(*@\label{cl_loadnew_begin}@*)
Ak{mu+1} = feval( Afun, kact+signed_mu );
Qk{mu+1} = eye(n);
fk{mu+1} = feval( ffun, kact+signed_mu );%(*@\label{cl_loadnew_end}@*)

for j=mu:-1:1 %(*@\label{cl_loop_begin}@*)
   % calculate a svd of E
   [U,S,V] = svd(Ek{1+j}); % Ek == U * S * V' %(*@\label{cl_svd_E_intheloop}@*)

   % determine the rank (only for sure)
   loc_r = sum( getdiag(S) > tol );
   if( loc_r ~= r(mu-j+1) )
      error(['r not invariant!  (k=',num2str(kact),',loc_r=',...
            num2str(loc_r),',r(',num2str(mu-j+1),')=',...
            num2str(r(mu-j+1)),')']);
   end

   % apply transformation
   S( (loc_r+1):m, : ) = zeros(m-loc_r,n);
   Ek{1+j} = S * V'; % = U'* Ekact
   Ak{1+j} = U'* Ak{1+j};
   fk{1+j} = U'* fk{1+j};

   % calculate a svd of A
   A2 = Ak{1+j}( loc_r+1:m, : );
   [U,S,V] = svd(A2); % Ak == U * S * V' %(*@\label{cl_svd_A_intheloop}@*)

   % determine the rank (only for sure)
   loc_h = sum( getdiag(S) > tol );
   if( loc_h ~= h(mu-j+1) )
      error(['h not invariant!  (k=',num2str(kact),',loc_h=',...
            num2str(loc_h),',h(',num2str(mu-j+1),')=',...
            num2str(h(mu-j+1)),')']);
   end

   % apply transformation
   S( loc_h+1:m-loc_r, : ) = zeros(m-loc_r-loc_h,n);
   Ek{j} = Ek{j}*V;
   Ak{1+j}( loc_r+1:m, : ) = S;
   Ak{1+j}( 1:loc_r, : )   = Ak{1+j}( 1:loc_r, : ) * V;
   fk{1+j}( loc_r+1:m, : ) = U'* fk{1+j}( loc_r+1:m, : );
   Qk{1+j} = Qk{1+j}*V;  %(*@\label{cl_remember_linmap}@*)

   % do block elimination (in A) %(*@\label{cl_elim_begin}@*)
   fk{1+j}(1:loc_r) = fk{1+j}(1:loc_r) - Ak{1+j}(1:loc_r,1:loc_h)*...
      (S(1:loc_h,1:loc_h)\fk{1+j}((loc_r+1):(loc_r+loc_h)));
   Ak{1+j}(1:loc_r,1:loc_h) = zeros(loc_r,loc_h);

   % do block elimination (in E)
   fk{j}(1:loc_r) = fk{j}(1:loc_r) + Ek{j}(1:loc_r,1:loc_h)*...
      (S(1:loc_h,1:loc_h)\fk{1+j}((loc_r+1):(loc_r+loc_h)));
   Ek{j}(1:loc_r,1:loc_h) = zeros(loc_r,loc_h); %(*@\label{cl_elim_end}@*)

end %(*@\label{cl_loop_end}@*)

[U,S,V] = svd(Ek{1}); % Ek == U * S * V' %(*@\label{cl_svd_E_outofloop}@*)

% determine the rank
loc_r = sum( getdiag(S) > tol );
if( length(r) >= mu+1 )
   if( loc_r ~= r(mu+1) )
      error(['r not invariant!  (k=',num2str(kact),',loc_r=',...
            num2str(loc_r),',r(',num2str(mu+1),')=',...
            num2str(r(mu+1)),')']);
   end
else
   r(mu+1) = loc_r;
end

% apply transformation
S( r(mu+1)+1:m, : ) = zeros(m-r(mu+1),n);
Ek{1} = S * V'; % = U'* Ekact
Ak{1} = U'* Ak{1};
fk{1} = U'* fk{1};

A2 = Ak{1}( r(mu+1)+1:m, : );
[U,S,V] = svd(A2); % Ak == U * S * V' %(*@\label{cl_svd_necces_outofloop}@*)

% determine the rank
loc_h = sum( getdiag(S) > tol );
if( length(h) >= mu+1 )
   if( loc_h ~= h(mu+1) )
      error(['h not invariant!  (k=',num2str(kact),',loc_h=',...
            num2str(loc_h),',h(',num2str(mu+1),')=',...
            num2str(h(mu+1)),')']);
   end
else
   h(mu+1) = loc_h;
end

% apply transformation
S( h(mu+1)+1:m-r(mu+1), : ) = zeros(m-r(mu+1)-h(mu+1),n);
Ek{1} = Ek{1};
Ak{1}( r(mu+1)+1:m, : ) = S;
Ak{1}( 1:r(mu+1), : ) = Ak{1}( 1:r(mu+1), : ) * V;
fk{1}( r(mu+1)+1:m, : ) = U'* fk{1}( r(mu+1)+1:m, : );
Qk{1} = Qk{1} * V;

% eliminate in A
fk{1}(1:r(mu+1)) = fk{1}(1:r(mu+1)) - Ak{1}(1:r(mu+1),1:h(mu+1))*...
   (S(1:h(mu+1),1:h(mu+1))\fk{1}((r(mu+1)+1):(r(mu+1)+h(mu+1))));
Ak{1}(1:r(mu+1),1:h(mu+1)) = zeros(r(mu+1),h(mu+1));
%(*@\label{cl_funct_end}@*)
% returns the diagonal elements of a matrix in a vector
function x=getdiag(A)

if( size(A,1) > 1 && size(A,2) > 1 )
   x = diag(A);
else
   if( prod(size(A)) > 0 )
      x = A(1,1);
   else
      x = [];
   end
end

function dumparray(Ek,Ak,fk,m,n) %(*@\label{cl_fct_dumparray}@*)
outsys = zeros(0,n);
outf   = [];
for i=1:prod(size(Ek))
   outsys = [ outsys , zeros((i-1)*m,n); ...
              zeros(m,(i-1)*n) , -Ak{i}, Ek{i} ];
   outf   = [outf;-fk{i}];
end
[outsys, outf]
