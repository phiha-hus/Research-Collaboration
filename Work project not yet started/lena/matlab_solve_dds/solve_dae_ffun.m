% Copyright (c) Tobias Bruell.  All rights reserved.
function f = loc_ffun(k)
global xxx_ffun;
global xxx_h;

f = feval(xxx_ffun, k*xxx_h, 0 );
