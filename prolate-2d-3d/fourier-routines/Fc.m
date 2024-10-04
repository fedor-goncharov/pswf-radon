% formard Fc via interpolated integral rule
% Fc : int_-1^1 exp(1i*c*x*y)f(y)dy

function val = Fc1d(f, f_grid, c, param, verbose=false, tol=1e-10, f_interp='linear')
   % Description:
   %   Fc in dimension 1 with interpolations and quadrature rules
   %   'f' is interpolated to a finer grid, so a care must be taken
   %   when representing 'f' on a grid
   % params:
   %   f     : values of function on a discrete grid
   %   f_grid : discrete grid inside [-1,1]
   %   f_interp : interpolation method (for interp1) for f
   %   c     : bandwidth parameter for PSWFs
   %   param : where to evaluate integal
   % output:
   %  scalar
   % assumed:
   %   f_grid is uniform
   % note:
   %  function calls many times interp1 (not efficient for large arrays)

   val = 0.0;
   A = f_grid(1); B=f_grid(end);
   if (A < -1.0 || B > 1.0)
     error('***** Fc1d -> domain of integrand exceeds [-1,1]. Abort.');
   endif

   % debug
   if verbose
     printf("%2.2f\n", (param - f_grid(1))/(f_grid(end) - f_grid(1)));
   endif

   exp_func = @(x) exp(1i*c*param*x);
   f_func = @(x) interp1(f_grid, f, x, f_interp);

   val = quadv(@(x) exp_func(x)*f_func(x), A, B, tol);
endfunction
