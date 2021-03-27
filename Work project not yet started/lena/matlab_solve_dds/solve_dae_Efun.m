% Copyright (c) Tobias Bruell.  All rights reserved.
function E = loc_Efun(k)
global xxx_Efun;
global xxx_h;
global xxx_discretization_method;

switch lower(xxx_discretization_method)
   case 'with_de'
      E = feval(xxx_Efun, (k+1)*xxx_h ,0 )/xxx_h;
   case 'ex_euler'
      E = feval(xxx_Efun, k*xxx_h ,0 )/xxx_h;
   otherwise
      error('Unknown discretization method');
end
