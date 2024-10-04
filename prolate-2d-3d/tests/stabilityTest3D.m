% testing PSWF inversion of integer order Hankel transforms
% agains additive Gaussian noise

amplitude = 1.0;
r1 = 0.25;
r2 = 0.5;
r3 = 0.75;
r4 = 1.0;
r5 = -1.0;
rend = -1.0;
phantom_func = @(x) amplitude*(x > r1 && x <= r2) ...
  + amplitude*(x > r3 && x <= r4) ...
  + amplitude*(x > r5 && x <= rend);

visualize_phantom1d = true;
if visualize_phantom1d
  figure;
  x = linspace(0, 1., 1000);
  y = arrayfun(phantom_func, x);
  plot(x,y)
  title("1d-phantom in experiment with Gaussian noise");
  xlabel ("x");
  ylabel ("y=f(x)");
endif% set main parameter 'c' for PSWFs

visualize_phantom2d = true;
if visualize_phantom2d
  figure;
  [XX,YY] = meshgrid(linspace(-1,1, 128), linspace(-1,1,128));
  RR = sqrt(XX.^2 + YY.^2);
  phantom_2d = arrayfun(@(r) phantom_func(r), RR);
  imshow(phantom_2d, [min(min(phantom_2d)), max(max(phantom_2d))]);
  title("phantom-2d in experiment with Gaussian noise");
endif

% set bandwidth
rsupp = 1.0; % support in [0, sigma] in the paper
bdwdth_restr = 10; % its restriction
bdwdth_nsize = 128; % 2048 - too much for 3D
bdwdth_grid = linspace(0, bdwdth_restr, bdwdth_nsize);
c = rsupp*bdwdth_restr;

% 1. NOISELESS RECONSTRUCTIONS -------------------------------------------------
%    find optimal number of modes

% generate [0, 5] orders up to 20 modes

l1_error_naive = zeros(6, 1);
l1_error_npswf = zeros(6, 21);
l2_error_naive = zeros(6, 1);
l2_error_npswf = zeros(6, 21);

recon_npswf_container = zeros(6, 21, 256);

common_grid = linspace(0.0, 1.0, 256);
phantom1d = arrayfun(phantom_func, linspace(0, 1, 256));
dx = common_grid(2) - common_grid(1);

do_replot_recon = true;
close_plots = true;
for n = 0:5
  printf("pre-selection job %d of %d\n", n, 5);
  fflush(stdout);

  order_hfinteger = n + 0.5;
  data_hfinteger = arrayfun(@(t) HankelIntegralTransformF(phantom_func, order_hfinteger, t, 0.0, rsupp), ...
    bdwdth_grid);

  % naive reconstruction via forward Hankel transform (it works for growing c)
  [recon_naive, recon_naive_grid] = reconFromHankelNaive(data_hfinteger, bdwdth_grid, c, order_hfinteger, 256, do_plotting = false, 128);

  l1_error_naive(n + 1) = trapz(abs(recon_naive - phantom1d))*dx;
  l2_error_naive(n + 1) = trapz((recon_naive - phantom1d).^2)*dx;

  % reconstruction using PSWFs
  for npswf=0:20

    [recon_npswf, recon_npswf_grid] = reconFromHankelPSWF3D(data_hfinteger, ...
                                            bdwdth_grid, ...
                                            c, ...
                                            uint8(n),
                                            uint8(npswf), ...
                                            256,
                                            256,
                                            256,
                                            "linear", ...
                                            0.0,
                                            0, ...
                                            false, ...
                                            128
                                           );

    recon_npswf_container(n + 1, npswf + 1, :) = recon_npswf;
    l1_error_npswf(n + 1, npswf + 1) = trapz(abs(recon_npswf - phantom1d))*dx;
    l2_error_npswf(n + 1, npswf + 1) = trapz((recon_npswf - phantom1d).^2)*dx;

    if (do_replot_recon == true)

      image_filename = sprintf("../prolate-2d-3d-experiment/images/3d-ord-%d-modes-%d.png", order_hfinteger, npswf);

      % phantom
      h = figure('position', [10 10 1762 256]);
      subplot(1, 4, 1);
      imagesc(phantom_2d, [min(min(phantom_2d)), max(max(phantom_2d))]);
      axis off;
      xlm = xlim;
      ylm = ylim;
      title("phantom", 'position', [mean(xlm) xlm(1)] ,'fontsize', 13);
      colorbar;
      set(gca, "fontsize", 12);

      % grid for reconstructions
      [YY, XX] = ndgrid(linspace(-1,1, 256), linspace(-1,1, 256));
      RR = sqrt(YY.^2 + XX.^2);

      % naive reconstruction
      subplot(1, 4, 2);
      recon_naive_3d = interp1(recon_naive_grid, recon_naive, RR, "linear", 0.0);
      imagesc(recon_naive_3d);
      axis off;
      colorbar;
      xlm = xlim;
      ylm = ylim;
      title(sprintf("3D naive reconstruction - order %2.1f", order_hfinteger), 'position', [mean(xlm) xlm(1)] , 'fontsize', 13);
      set(gca, "fontsize", 12);

      % npswf reconstruction
      subplot(1, 4, 3);
      recon_npswf_3d = interp1(recon_npswf_grid, recon_npswf, RR, "linear", 0.0);
      imagesc(real(recon_npswf_3d));
      axis off;
      colorbar;
      xlm = xlim;
      ylm = ylim;
      title(sprintf("3D reconstruction - order %2.1f, %d modes", order_hfinteger, npswf), 'position', [mean(xlm) xlm(1)], 'fontsize', 13);
      set(gca, "fontsize", 12);

      % central slice
      subplot(1, 4, 4);
      plot(recon_npswf_grid, real(recon_npswf), sprintf(";noiseless PSWF;", npswf), "linewidth", 1.5);
      hold on;
      plot(recon_npswf_grid, arrayfun(phantom_func, recon_npswf_grid), ";phantom;", "linewidth", 1.5);
      hold on;
      plot(recon_naive_grid, recon_naive, ";naive;", "linewidth", 1.5);
      hold off;
      xlm = xlim;
      ylm = ylim;
      title(sprintf("1D central slice - order %2.1f, %d modes", order_hfinteger, npswf), 'position', [mean(xlm) ylm(2)], 'fontsize', 13);
      legend("location", "northeastoutside");
      set(gca, "fontsize", 12);

      print(image_filename, "-dpng", "-color", "-S1762,256");
      %saveas(h, image_filename);
      if (close_plots == true)
        close(h);
      endif
    endif
  endfor
endfor

save stability-test-3d.mat recon_npswf_container l1_error_naive l1_error_npswf l2_error_naive l2_error_npswf



##% ----------------------------------------------------------------------------
##
##
% 2. ----------------------- l1/l2-errors plots --------------------------------
do_plot_errors = false;
for n = 0:5

  ord = n + 0.5;
  printf("Order %d\n", i);
  [l1_min_vals l1_min_inds] = min(l1_error_npswf(n+1, :)(:));
  [l2_min_vals l2_min_inds] = min(l2_error_npswf(n+1, :)(:));
  printf("\tL1-argmin at first mode : %d\n", l1_min_inds(1)-1);
  printf("\tL2-argmin at first mode : %d\n", l2_min_inds(1)-1);

  if (do_plot_errors == true)
    % l1
    h=figure('position', [10 10 722 522]);
    image_filename = sprintf("../prolate-2d-3d-experiment/images/errors/3d-l1-errs-ord-%2.1f.png", ord);

    plot(0:20, l1_error_naive(n+1)*ones(21, 1), ";naive;", "linewidth", 1.5);
    title(sprintf("L1-errors for order %f", ord), 'fontsize', 13);
    xlim([0 20]);
    ylim([0 1]);
    hold on;
    plot(0:20, l1_error_npswf(n+1, :), ";2D npswf reconstructions;", "linewidth", 1.5);
    legend("location", "northeastoutside");
    set(gca, "fontsize", 12);

    saveas(h, image_filename);
    close(h);

    % l2
    h=figure('position', [10 10 722 522]);
    image_filename = sprintf("../prolate-2d-3d-experiment/images/errors/3d-l2-errs-ord-%2.1f.png", ord);

    plot(0:20, l2_error_naive(n+1)*ones(21, 1), ";naive;", "linewidth", 1.5);
    title(sprintf("L2-errors for order %2.1f", ord), 'fontsize', 13);
    xlim([0 20]);
    ylim([0 1]);
    hold on;
    plot(0:20, l2_error_npswf(n+1, :), ";2D npswf reconstructions;", "linewidth", 1.5);
    legend("location", "northeastoutside");
    set(gca, "fontsize", 12);

    saveas(h, image_filename);
    close(h);
  endif
endfor

% ------------------------------------------------------------------------------


##% 3. NOISE STABILITY TEST FOR MINIMAL ORDERS ---------------------------------
##
##noise_levels = [0.01 0.05 0.07 0.1 0.12 0.15];
##order_integer = 4.0;
##data_integer = arrayfun(@(t) HankelIntegralTransformF(phantom_func, order_integer, t, 0.0, rsupp), ...
##    bdwdth_grid);
##
##dc = bdwdth_grid(2) - bdwdth_grid(1);
##l2_data_norm = sqrt(trapz(data_integer.^2)*dc);
##npswf = 12; % minimal number of modes for order 0
##
##for r=noise_levels
##  r
##  % setting noise strength for predefined relative L2-error level
##  log_sigma = log(r)+log(l2_data_norm)-log(sqrt(2))-(lgamma(129.0/2.0)-lgamma(64))-0.5*log(dc);
##  sigma = exp(log_sigma);
##
##  % check that noise mapping is correct
##  % z = sigma*randn(size(bdwdth_grid));
##  % sqrt(dc*trapz(z.^2))/l2_data_norm
##
##  nsamples = 1000;
##  z = sigma*randn(nsamples, 128);
##  recon_noisy = zeros(nsamples, 256);
##
##  [recon_npswf_noiseless, recon_npswf_grid_noiseless] = reconFromHankelPSWF2D(hankel_data=data_integer, ...
##                                                             hankel_grid=bdwdth_grid, ...
##                                                             c=c, ...
##                                                             order=uint8(order_integer),
##                                                             npswf=uint8(npswf), ...
##                                                             phi_grid_size=256,
##                                                             recon_grid_size=256,
##                                                             interp_method = "linear", ...
##                                                             rpadd = 1.3,
##                                                             verbose_level = 0, ...
##                                                             do_plotting = false, ...
##                                                             plot2d_size = 64
##                                                            );
##
##  for i = 1:nsamples
##    i
##    data_integer_noisy = data_integer + z(i, :);
##    [recon_npswf, recon_npswf_grid] = reconFromHankelPSWF2D(hankel_data=data_integer_noisy, ...
##                                                             hankel_grid=bdwdth_grid, ...
##                                                             c=c, ...
##                                                             order=uint8(order_integer),
##                                                             npswf=uint8(npswf), ...
##                                                             phi_grid_size=256,
##                                                             recon_grid_size=256,
##                                                             interp_method = "linear", ...
##                                                             rpadd = 1.3,
##                                                             verbose_level = 0, ...
##                                                             do_plotting = false, ...
##                                                             plot2d_size = 64
##                                                            );
##    recon_noisy(i, :) = recon_npswf(:);
##  endfor
##
##  recon_mean = mean(recon_noisy, 1);
##  recon_std = std(recon_noisy, 1);
##
##  h=figure('position', [10 10 822 422]);
##  image_filename = sprintf("../prolate-2d-3d-experiment/images/errors/errs-bands-ord-%d-modes-%d-noise-%2.2f.png", order_integer, npswf, r*100);
##
##  plot(common_grid, recon_mean, sprintf(";mean, noiselevel %2.1f%%;", 100*r), "linewidth", 1.5);
##  ylm = ylim;
##  xlm = xlim;
##  title(sprintf("95%% error bands for npswf, order=%d, modes=%d", order_integer, npswf), 'position', [mean(xlm) ylm(2)]);
##  legend("autoupdate", "off");
##  patch([common_grid flip(common_grid)], [recon_mean-2*recon_std flip(recon_mean + 2*recon_std)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
##  legend("autoupdate", "on");
##  ylim(ylm);
##  hold on;
##  plot(recon_npswf_grid_noiseless, recon_npswf_noiseless, ";npswf noiseless;", "linewidth", 1.5);
##  %hold on;
##  plot(recon_npswf_grid, arrayfun(phantom_func, recon_npswf_grid), ";phantom;", "linewidth", 1.5);
##
##  legend("location", "northeastoutside");
##  set(gca, 'fontsize', 13);
##  saveas(h, image_filename);
##endfor
##
##%-------------------------------------------------------------------------------



