% reconstruct function from sinogram data in 2D
% standard 'iradon' (filtered backprojection) on uniform grid

function [v_data, v_grid_x, v_grid_y] = TransformSinogramToV2D(sinogram_data, ...
                                                               sinogram_s_grid, ...
                                                               sinogram_phi_grid,
                                                               padd_size=128)
  % Description:
  %   Radon inversion in 2D  from 'sinogram_data'.
  %   Can accept complex values.
  % params:
  %
  %  padd_size : padd symmetrically with 'padd_size' zeros for each projection
  %   ...
  % assumed:
  %   phi spans [0, pi], endpoints included

  % padd with zeros for stability
  sinogram_padded = padarray(sinogram_data, [0, padd_size]);
  scale = ceil(size(sinogram_data)(2)/2); % for 'iradon'

  % scale correctly - needed for later 'iradon'
  sinogram_padded_real = real(sinogram_padded)*scale;
  sinogram_padded_imag = imag(sinogram_padded)*scale;

  % transpose sinograms and feed to standard 'iradon'
  theta = linspace(0, 1, size(sinogram_data)(1))*180; % DEBUG: sinogram_phi_grid not used, remove endpoint
  recon_real = iradon(sinogram_padded_real', theta);
  recon_imag = iradon(sinogram_padded_imag', theta);

  v_data = recon_real + 1j*recon_imag;

  extended_limit = size(sinogram_padded)(2)/(size(sinogram_data)(2)*sqrt(2));
  [v_grid_x,v_grid_y] = meshgrid(linspace(-extended_limit, extended_limit, size(v_data)(1)),
                   linspace(-extended_limit, extended_limit, size(v_data)(1)));

endfunction
