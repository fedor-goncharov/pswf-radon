% functional Hankel integral transform
% accepts calling function

function val = HankelIntegralTransformF(f, ord, arg, supp_fa, supp_fb, tol=1e-10)
  % no checks of user input

  % f     : function
  % ord   : order of Bessel function
  % arg   : transform evaluated at point
  % supp_fa, supp_fb : support interval of 'f'

  integrand = @(x) f(x)*besselj(ord, arg*x)*sqrt(arg*x);
  val = quadv(integrand, supp_fa, supp_fb, tol);
endfunction


