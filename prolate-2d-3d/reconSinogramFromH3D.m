% reconstruct 3D sinogram from 1D H-data specifically
% for 1D embedded data
%   this equivalent to combination of 'TransformHToFcSinogram3D' and 'TransformFcSinogramToSinogram3D'
%   but works faster is it does not need to repeat same (up to multiplicative factor)
%   reconstructions

function [sinogram_data, ...
          sinogram_theta_grid, ...
          sinogram_phi_grid, ...
          sinogram_s_grid] = reconSinogramFromH3D(h_data, ...
                                                  h_grid, ...
                                                  theta_grid_size, ...
                                                  phi_grid_size, ...
                                                  c, ...
                                                  npswf, ...
                                                  order, ...
                                                  rsupp, ...
                                                  verbose_level=1
                                                  )
  % Fc-inversion
  if (verbose_level >= 1)
    printf("reconSinogramFromH3D :: inversion of 1D signal using PSWFs with %d modes...", npswf);
  endif

  verbose_invFc = false;
  if (verbose_level >= 2)
    verbose_invFc = true;
  endif
  sinogram_const_slice = invFc(h_data, h_grid, c, npswf, verbose = verbose_invFc);
  if (verbose_level >= 1)
    printf("Done.\n");
  endif

  % make 3D sinogram from slice
  if (verbose_level >= 1)
    printf("reconSinogramFromH3D :: embedding in a 3D sinogram...");
  endif
  n = order-0.5;
  phi_grid = linspace(0, 2*pi, phi_grid_size + 1)(1:end-1); % no endpoint due to periodicity
  theta_grid = acos(linspace(-1, 1, theta_grid_size + 2)(2:end-1)); % exclude 'north' and 'south' poles

  [PP, TT, SS] = meshgrid(phi_grid, theta_grid, h_grid); % 2D grid of angles (theta, phi)
  sinogram_phi_grid = PP;
  sinogram_theta_grid = TT;
  sinogram_s_grid = SS;

  scalar_coeff = (2*pi)^(1.5)*(1j)^n/(rsupp^3); % set scalar multiplier
  sinogram_const_slice_rep = squeeze(permute(repmat(sinogram_const_slice, [1 1 theta_grid_size phi_grid_size]), [1 3 4 2])); % get elementwise product 'h'x'harmonicY'
  sinogram_data = scalar_coeff*harmonicY(n, 0, TT, PP, 'norm', true) .* sinogram_const_slice_rep;

  if (verbose_level >= 1)
    printf("Done.\n");
  endif

endfunction
