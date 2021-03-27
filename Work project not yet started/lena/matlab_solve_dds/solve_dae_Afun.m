% Copyright (c) Tobias Bruell.  All rights reserved.
function A = loc_Afun(k)
global xxx_Afun;
global xxx_Efun;
global xxx_h;
global xxx_discretization_method;

switch lower(xxx_discretization_method)
   case 'with_de'
      A = feval(xxx_Afun, k*xxx_h, 0 ) + ...
          feval(xxx_Efun, k*xxx_h, 0 )/xxx_h + ...
          feval(xxx_Efun, k*xxx_h, 1 );
   case 'ex_euler'
      A = feval(xxx_Afun, k*xxx_h, 0 ) + ...
          feval(xxx_Efun, k*xxx_h, 0 )/xxx_h;
   otherwise
      error('Unknown discretization method');
end
