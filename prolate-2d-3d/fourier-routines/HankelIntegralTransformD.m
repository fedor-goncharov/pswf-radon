% Hankel integral transform via discrete integration

function val = HankelIntegralTransformD(f_data, f_grid, ord, arg)
  % no checks of user input

  % f_data     : array of values
  % f_grid     : grid of integration
  % ord   : order of Bessel function
  % arg   : transform evaluated at point
  % supp_fa, supp_fb : support interval of 'f'

  integrand = f_data.*besselj(ord, arg*f_grid).*sqrt(arg*f_grid);
  dx = f_grid(2) - f_grid(1);
  val = trapz(integrand)*dx;
  %val = IntegrateSimpson1D(integrand, f_grid(2)-f_grid(1));
endfunction


