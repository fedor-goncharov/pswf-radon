% transform FcSinogram to Sinogram
% see formulas in Theorem 2.1 in the original text

function [sinogram_data, sinogram_s_grid, sinogram_phi_grid] = TransformFcSinogramToSinogram2D(fcsinogram_data, ...
                                                                        fcsinogram_h_grid, ...
                                                                        fcsinogram_phi_grid, ...
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
  sinogram_s_grid = fcsinogram_h_grid;
  sinogram_phi_grid = fcsinogram_phi_grid;

  % loop for phi and invFc on each
  if verbose
    printf("TransformFcSinogramToSinogram -> Size of sinogram: %d x %d\n", size(fcsinogram_data)(1), size(fcsinogram_data)(2));
  endif
  for i=1:size(fcsinogram_data)(1)
    if verbose
      printf("TransformFcSinogramToSinogram -> Projection %d\n", i);
    endif
    sinogram_data(i, :) = invFc(fcsinogram_data(i, :), fcsinogram_h_grid(i, :), c, npswf, verbose);
  endfor

endfunction
