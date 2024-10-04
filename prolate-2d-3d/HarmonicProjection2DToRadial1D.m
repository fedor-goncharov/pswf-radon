% average 2D complex-valued signal over circular harmonics on S^1
% computes integrals f(r, phi)exp(ik\phi) d\phi on the circle for various r
% via standard trapezoidal rule


function [rad_proj, rad_proj_grid] = HarmonicProjection2DToRadial1D(f_data, ...
                                                              f_grid_x, ...
                                                              f_grid_y, ...
                                                              rad_proj_grid_size, ...
                                                              order, ...
                                                              phi_grid_size, ...
                                                              interp_method = "linear", ...
                                                              divide_2pi = true, ...
                                                              grid_pos_o = true
                                                              )
  rad_proj = zeros(1, rad_proj_grid_size);
  rad_proj_grid = linspace(0, 1, rad_proj_grid_size);
  phi_grid = linspace(0, 2*pi, phi_grid_size);

  harmonic_vals = exp(-1j*order*phi_grid);
  phi_o_factor = 1.; % orientation factor
  if (grid_pos_o == false)
    phi_o_factor = -1.;
  endif

  for r_ind = 1:rad_proj_grid_size
    r = rad_proj_grid(r_ind);
    % DEBUG - changed for NFFT
    xi = r*cos(phi_o_factor*phi_grid); % note '-phi' - this is because grids 'x','y' are [-1,1] are in inverse order
    yi = r*sin(phi_o_factor*phi_grid);
    f_r_vals = interp2(f_grid_x, f_grid_y, f_data, xi, yi, interp_method, 0.0);
    rad_proj(r_ind) = trapz(f_r_vals.*harmonic_vals)*(phi_grid(2)-phi_grid(1));
  endfor
  if (divide_2pi == true)
    rad_proj /= 2*pi;
  endif
endfunction
