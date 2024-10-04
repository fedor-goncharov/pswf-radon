function vals = TransformFcSimpson1d(f, f_grid, c, params, verbose=false)
   % Description:
   %   Fc in dimension 1 with Simpson integration rule
   % params:
   %   f     : values of function on a discrete grid
   %   f_grid : discrete grid inside [-1,1]
   %   f_interp : interpolation method (for interp1) for f
   %   c     : bandwidth parameter for PSWFs
   %   param : array where to evaluate integal
   % output:
   %  scalar
   % assumed:
   %   f_grid is uniform
   % note:
   %  function calls many times interp1 (not efficient for large arrays)

   vals = zeros(1, length(params));
   if (f_grid(1) < -1.0 || f_grid(end) > 1.0)
     error('***** FcSimpson1d -> domain of integrand exceeds [-1,1]. Abort.');
   endif

   dx = f_grid(2) - f_grid(1);
   for ind=1:length(params)
     if verbose
       ind
     endif
     param = params(ind);
     k = arrayfun(@(x) exp(1i*c*param*x), f_grid);
     vals(ind) = IntegrateSimpson1D(f .* k, dx);
   endfor
endfunction
