% transform H-data to Fc transform of the sinogram in 3D
% see formulas (2.7) in the original article

function [fcsinogram_data, ...
          fcsinogram_theta_grid, ...
          fcsinogram_phi_grid, ...
          fcsinogram_s_grid] = TransformHToFcSinogram3D(h_data, h_grid, theta_grid_size, phi_grid_size, order, rsupp)
  % params:
  %   h_data   : h-transform data
  %   h_grid   : h-transform grid
  %   phi_grid : array grid of phi in radians (assumed not vec)
  %   rsupp    : radius of signal support
  % return:
  %   1 sinogram (complex), 3 grids (real; h-grid, phi-grid, theta-grid) all
  %   from meshgrid with order (theta, phi, h)
  %   grid in theta is not uniform but 'acos' one
  % dependencies:
  %   harmonicY
  % note:
  %   order is half-integer
  %   h_grid uniform in [-1,1] endpoints included
  %   phi_grid uniform in [0, 2*pi] endpoints included

  % set grids
  n = order-0.5;
  phi_grid = linspace(0, 2*pi, phi_grid_size);
  theta_grid = acos(linspace(-1, 1, theta_grid_size + 2)(2:end-1)); % exclude zero and pi

  [PP, TT, SS] = meshgrid(phi_grid, theta_grid, h_grid); % 2D grid of angles (theta, phi)
  fcsinogram_phi_grid = PP;
  fcsinogram_theta_grid = TT;
  fcsinogram_s_grid = SS;
  % set scalar multiplier
  scalar_coeff = (2*pi)^(1.5)*(1j)^n/(rsupp^3);
  % get elementwise product 'h'x'harmonicY'
  h_data_rep = squeeze(permute(repmat(h_data, [1 1 length(theta_grid) length(phi_grid)]), [1 3 4 2]));
  fcsinogram_data = scalar_coeff*harmonicY(n, 0, TT, PP, 'norm', true) .* h_data_rep;
endfunction
