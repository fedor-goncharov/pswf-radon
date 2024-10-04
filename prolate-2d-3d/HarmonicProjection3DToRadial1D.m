% average 3D complex-valued signal over circular harmonics on S^2 for various r
% via standard trapezoidal rule with linear interpolation over 3D grid


function [rad_proj, rad_proj_grid] = HarmonicProjection3DToRadial1D(f_data, ...
                                                              f_grid_x, ...
                                                              f_grid_y, ...
                                                              f_grid_z,
                                                              rad_proj_grid_size, ...
                                                              n, ...
                                                              phi_grid_size, ...
                                                              theta_grid_size, ...
                                                              interp_method = "linear", ...
                                                              from_nfft = false,
                                                              verbose_level = 1
                                                              )
  % assumed:
  %   'n' is integer
  %   f_grid_x, .. f_grid_z are ordered monotonically, derived from 'meshgrid' or 'ngrid'

  % reconstruction grid
  rad_proj = zeros(1, rad_proj_grid_size);
  rad_proj_grid = linspace(0, 1, rad_proj_grid_size);

  % set spherical harmonics
  phi_grid = linspace(0, 2*pi, phi_grid_size + 1)(1:end-1);
  dphi = phi_grid(2) - phi_grid(1);
  u_grid = linspace(-1., 1., theta_grid_size + 2)(2:end-1);
  du = u_grid(2) - u_grid(1);
  theta_grid = acos(u_grid);
  [PP, TT] = meshgrid(phi_grid, theta_grid);
  harmonic_vals = harmonicY(n, 0, TT, PP, 'norm', true)(:);

  % data reconstructed from NUFFT must be reordered (original order [z, y, x])
  if from_nfft == true
    f_data = permute(f_data, [3 2 1]);
  endif
  for r_ind = 1:rad_proj_grid_size
    if verbose_level >= 2
      printf("\titeration %d of (%d)\n", r_ind, rad_proj_grid_size);
    endif

    r = rad_proj_grid(r_ind);
    xi = (r*cos(PP).*sin(TT))(:);
    yi = (r*sin(PP).*sin(TT))(:);
    zi = (r*cos(TT))(:);

    f_r_vals = interpn(f_grid_x, f_grid_y, f_grid_z, f_data, xi, yi, zi, interp_method, 0.0);
    rad_proj(r_ind) = sum(f_r_vals.*harmonic_vals)*dphi*du;
  endfor
endfunction
