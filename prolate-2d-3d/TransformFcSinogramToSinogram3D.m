% transform FcSinogram to Sinogram
% see formulas in Theorem 2.1 in the original text

function [sinogram_data, ...
          sinogram_theta_grid, ...
          sinogram_phi_grid, ...
          sinogram_s_grid] = TransformFcSinogramToSinogram3D( ...
                                                              fcsinogram_data, ...
                                                              fcsinogram_theta_grid, ...
                                                              fcsinogram_phi_grid, ...
                                                              fcsinogram_s_grid, ...
                                                              c, ...
                                                              npswf,...
                                                              verbose = false
                                                              )
  % Description
  %   inversion of Fc for each projection in sinogram data
  % params:
  %   fcsinogram_data : sinogram 2d-matrix
  %   fcsinogram_h_grid : scalar shifts matrix (from meshgrid)
  %   fcsingoram_phi_grid : projection angles matrix (from meshgrid)
  %   c : bandwidth parameter for PSWFs
  %   npswf :  number of PSWFs used
  %   verbose : print additional info
  % output:
  %   inverse Fc for each projection
  % assumed:
  %   sinogram has size length(phi_grid) x length(h_grid)

  % init containers
  sinogram_data = zeros(size(fcsinogram_data));
  sinogram_theta_grid = fcsinogram_theta_grid;
  sinogram_s_grid = fcsinogram_s_grid;
  sinogram_phi_grid = fcsinogram_phi_grid;

  % loop for phi and theta invFc on each
  if verbose
    printf("TransformFcSinogramToSinogram -> Size of sinogram: %d x %d x %d\n", size(fcsinogram_data)(1), size(fcsinogram_data)(2), size(fcsinogram_data)(3));
  endif
  for ind_theta=1:size(fcsinogram_data)(1)
    if verbose
      printf("TransformFcSinogramToSinogram -> Projection %d\n", ind_theta);
    endif
    for ind_phi=1:size(fcsinogram_data)(2)
      % reshape inputs to row vectors
      fcsinogram_data_slice = reshape(fcsinogram_data(ind_theta, ind_phi, :), [1 size(fcsinogram_data)(end)]);
      fcsinogram_s_grid_slice = reshape(fcsinogram_s_grid(ind_theta, ind_phi,:), [1 size(fcsinogram_s_grid)(end)]);
      sinogram_data(ind_theta, ind_phi, :) = invFc( ...
                                                   fcsinogram_data_slice, ...
                                                   fcsinogram_s_grid_slice, ...
                                                   c, npswf, verbose ...
                                                   );
    endfor
  endfor
endfunction
