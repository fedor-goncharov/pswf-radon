% reconstruct 3D sinogram from 1D H-data specifically
% for 1D embedded data
%   this equivalent to combination 'TransformHToFcSinogram2D' and 'TransformFcSinogramToSinogram2D'
%   but works faster is it does not need to repeat same (up to multiplicative factor)
%   reconstructions

function [sinogram_data, ...
          sinogram_phi_grid, ...
          sinogram_s_grid] = reconSinogramFromH2D(h_data, ...
                                                  h_grid, ...
                                                  phi_grid_size, ...
                                                  c, ...
                                                  npswf, ...
                                                  order, ...
                                                  rsupp, ...
                                                  verbose_level=1
                                                  )
  % Fc-inversion
  if (verbose_level >= 1)
    printf("reconSinogramFromH2D :: inversion of 1D signal using PSWFs with %d modes...", npswf);
  endif

  verbose_invFc = false;
  if (verbose_level >= 2)
    verbose_invFc = true;
  endif
  sinogram_const_slice = invFc(h_data, h_grid, c, npswf, verbose = verbose_invFc);
  if (verbose_level >= 1)
    printf("Done.\n");
  endif

  % make 2D sinogram from slice
  if (verbose_level >= 1)
    printf("reconSinogramFromH3D :: embedding in a 3D sinogram...");
  endif
  n = order;
  phi_grid = linspace(0, 2*pi, phi_grid_size + 1)(1:end-1); % no endpoint due to periodicity

  [PP, SS] = ndgrid(phi_grid, h_grid);
  sinogram_phi_grid = PP;
  sinogram_s_grid = SS;

  scalar_coeff = (2*pi)*(1j)^n/(rsupp^2); % set scalar multiplier
  sinogram_const_slice_rep = squeeze(permute(repmat(sinogram_const_slice, [1 1 phi_grid_size]), [1 3 2]));
  sinogram_data = scalar_coeff*exp(1j*n*PP) .* sinogram_const_slice_rep;

  if (verbose_level >= 1)
    printf("Done.\n");
  endif

endfunction
