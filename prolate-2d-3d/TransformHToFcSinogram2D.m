% transform H-data to Fc transform of the sinogram in 2D
% see formula (2.5) in the original article

function [fcsinogram_data, fcsinogram_h_grid, fcsinogram_phi_grid] = TransformHToFcSinogram2D(h_data, ...
                                                                                              h_grid, ...
                                                                                              phi_grid_size, ...
                                                                                              order, ...
                                                                                              rsupp, ...
                                                                                              to2pi = false)
  % params:
  %   h_data   : h-transform data
  %   h_grid   : h-transform grid
  %   phi_grid : array grid of phi in radians (assumed not vec)
  %   rsupp    : radius of signal support
  % Note:
  %   order is assumed integer
  %   h_grid uniform in [-1,1] endpoints included
  %   phi_grid uniform in [0, pi] endpoints included


  % set 'phi' grid
  if (to2pi == true)
    phi_grid = linspace(0, 2*pi, phi_grid_size + 1)(1:phi_grid_size);
  else
    phi_grid = linspace(0, pi, phi_grid_size + 1)(1:phi_grid_size);
  endif

  % set data
  coeff = 2*pi*(1j)^order/(rsupp^2);
  fcsinogram_data = coeff*vec(exp(1j*order*phi_grid))*h_data;

  % set complete grid
  [fcsinogram_h_grid, fcsinogram_phi_grid] = meshgrid(h_grid, phi_grid);
endfunction
