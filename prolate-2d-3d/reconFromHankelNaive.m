% reconstruction from Hankel transform limited on [0, c]
% via second Hankel transform (with zero extension)

function [naive_recon, recon_grid] = reconFromHankelNaive(hankel_data, ...
                                                   hankel_grid, ...
                                                   c, ...
                                                   order, ...
                                                   recon_grid_size, ...
                                                   do_plotting = false, ...
                                                   plot2d_size = 64
                                                  )
  % Description:
  %   ...
  % params:
  %   ...
  % assumed:
  %   reconstructed signal is supported on [0, 1]

  assert(hankel_grid(1) == 0.0 && hankel_grid(end) == c);

  % take Hankel transform through integration on provided grid (no interpolation)
  recon_grid = linspace(0, 1, recon_grid_size);
  naive_recon = arrayfun(@(t) HankelIntegralTransformD(hankel_data, hankel_grid, order, t), recon_grid);
  if (do_plotting == true)
    figure, plot(linspace(0, 1, recon_grid_size), naive_recon), title(sprintf("1D naive reconstruction via Hankel transform, order=%.2f, c=%.2f", order, c));
  endif
  [XX, YY] = meshgrid(linspace(-1,1, plot2d_size), linspace(-1,1, plot2d_size));
  RR = sqrt(XX.^2 + YY.^2);
  if (do_plotting == true)
    naive_recon_2d = interp1(recon_grid, naive_recon, RR, "nearest", 0.0);
    figure, imagesc(naive_recon_2d, [min(min(naive_recon_2d)), max(max(naive_recon_2d))]), ...
      title(sprintf("2D naive reconstruction via Hankel transform, order=%.2f, c=%.2f", order, c)), colorbar;
  endif
endfunction

