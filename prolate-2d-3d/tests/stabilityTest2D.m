% testing PSWF inversion of integer order Hankel transforms
% agains additive Gaussian noise

amplitude = 1.0;
r1 = 0.25;
r2 = 0.5;
r3 = 0.75;
r4 = 0.95;
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
bdwdth_nsize = 128; % 2048
bdwdth_grid = linspace(0, bdwdth_restr, bdwdth_nsize);
c = rsupp*bdwdth_restr;

% 1. NOISELESS RECONSTRUCTIONS -------------------------------------------------
%    Comparison phantom vs. naive vs. npswf vs. 1d-section

% generate [0, 5] orders up to 20 modes

l1_error_naive = zeros(6, 1);
l1_error_npswf = zeros(6, 21);
l2_error_naive = zeros(6, 1);
l2_error_npswf = zeros(6, 21);

common_grid = linspace(0.0, 1.0, 256);
phantom1d = arrayfun(phantom_func, linspace(0, 1, 256));
dx = common_grid(2) - common_grid(1);

do_replot_recon = true;
for order_integer = 4:4

  data_integer = arrayfun(@(t) HankelIntegralTransformF(phantom_func, order_integer, t, 0.0, rsupp), ...
    bdwdth_grid);

  % naive reconstruction via forward Hankel transform (it works for growing c)
  [recon_naive, recon_naive_grid] = reconFromHankelNaive(data_integer, bdwdth_grid, c, order_integer, 256, do_plotting = false, 128);

  l1_error_naive(order_integer + 1) = trapz(abs(recon_naive - phantom1d))*dx;
  l2_error_naive(order_integer + 1) = trapz((recon_naive - phantom1d).^2)*dx;

  % reconstruction using PSWFs
  for npswf=4:4
    [recon_npswf, recon_npswf_grid] = reconFromHankelPSWF2D(hankel_data=data_integer, ...
                                                             hankel_grid=bdwdth_grid, ...
                                                             c=c, ...
                                                             order=uint8(order_integer),
                                                             npswf=uint8(npswf), ...
                                                             phi_grid_size=256,
                                                             recon_grid_size=256,
                                                             interp_method = "linear", ...
                                                             rpadd = 2.0,
                                                             verbose_level = 1, ...
                                                             do_plotting = false, ...
                                                             plot2d_size = 64
                                                            );
    l1_error_npswf(order_integer + 1, npswf + 1) = trapz(abs(recon_npswf - phantom1d))*dx;
    l2_error_npswf(order_integer + 1, npswf + 1) = trapz((recon_npswf - phantom1d).^2)*dx;

    if (do_replot_recon == true)

      image_filename = sprintf("../prolate-2d-3d-experiment/images/2d-ord-%d-modes-%d.png", order_integer, npswf);

      % grid
      [YY, XX] = ndgrid(linspace(-1,1, 256), linspace(-1,1, 256));
      RR = sqrt(YY.^2 + XX.^2);

      % phantom
      h = figure('position', [10 10 1762 292]);
      subplot(1, 4, 1);
      imagesc(phantom_2d, [min(min(phantom_2d)), max(max(phantom_2d))]);
      axis off;
      xlm = xlim;
      ylm = ylim;
      title("phantom", 'position', [mean(xlm) xlm(1)] ,'fontsize', 13);
      colorbar;
      set(gca, "fontsize", 12);

      % naive reconstruction
      subplot(1, 4, 2);
      recon_naive_2d = interp1(recon_naive_grid, recon_naive, RR, "linear", 0.0);
      imagesc(recon_naive_2d);
      axis off;
      colorbar;
      xlm = xlim;
      ylm = ylim;
      title(sprintf("2D naive reconstruction - order %d", order_integer), 'position', [mean(xlm) xlm(1)] , 'fontsize', 13);
      set(gca, "fontsize", 12);

      % npswf reconstruction
      subplot(1, 4, 3);
      recon_npswf_2d = interp1(recon_npswf_grid, recon_npswf, RR, "linear", 0.0);
      imagesc(real(recon_npswf_2d));
      axis off;
      colorbar;
      xlm = xlim;
      ylm = ylim;
      title(sprintf("2D reconstruction - order %d, %d modes", order_integer, npswf), 'position', [mean(xlm) xlm(1)], 'fontsize', 13);
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
      title(sprintf("1D central slice - order %d, %d modes", order_integer, npswf), 'position', [mean(xlm) ylm(2)], 'fontsize', 13);
      legend("location", "northeastoutside");
      set(gca, "fontsize", 12);

      %saveas(h, image_filename);
      %close(h);
    endif
  endfor
endfor

% ------------------------------------------------------------------------------


% 2. l1/l2-errors plots ---------------------------------------------------------
do_plot_errors = false;
for i = 0:5


  printf("Order %d\n", i);
  [l1_min_vals l1_min_inds] = min(l1_error_npswf(i+1, :)(:));
  [l2_min_vals l2_min_inds] = min(l2_error_npswf(i+1, :)(:));
  printf("\tL1-argmin at first mode : %d\n", l1_min_inds(1)-1);
  printf("\tL2-argmin at first mode : %d\n", l2_min_inds(1)-1);

  if (do_plot_errors == true)
    % l1
    h=figure('position', [10 10 722 522]);
    image_filename = sprintf("../prolate-2d-3d-experiment/images/errors/l1-errs-ord-%d.png", i);

    plot(0:20, l1_error_naive(i+1)*ones(21, 1), ";naive;", "linewidth", 1.5);
    title(sprintf("L1-errors for order %d", i), 'fontsize', 13);
    xlim([0 20]);
    ylim([0 1]);
    hold on;
    plot(0:20, l1_error_npswf(i+1, :), ";2D npswf reconstructions;", "linewidth", 1.5);
    legend("location", "northeastoutside");
    set(gca, "fontsize", 12);

    saveas(h, image_filename);
    close(h);

    % l2
    h=figure('position', [10 10 722 522]);
    image_filename = sprintf("../prolate-2d-3d-experiment/images/errors/l2-errs-ord-%d.png", i);

    plot(0:20, l2_error_naive(i+1)*ones(21, 1), ";naive;", "linewidth", 1.5);
    title(sprintf("L2-errors for order %d", i), 'fontsize', 13);
    xlim([0 20]);
    ylim([0 1]);
    hold on;
    plot(0:20, l2_error_npswf(i+1, :), ";2D npswf reconstructions;", "linewidth", 1.5);
    legend("location", "northeastoutside");
    set(gca, "fontsize", 12);

    saveas(h, image_filename);
    close(h);
  endif
endfor

% ------------------------------------------------------------------------------


% 3. NOISE STABILITY TEST FOR MINIMAL ORDERS -----------------------------------

noise_levels = [0.01 0.05 0.07 0.1 0.12 0.15];
order_integer = 4.0;
data_integer = arrayfun(@(t) HankelIntegralTransformF(phantom_func, order_integer, t, 0.0, rsupp), ...
    bdwdth_grid);

dc = bdwdth_grid(2) - bdwdth_grid(1);
l2_data_norm = sqrt(trapz(data_integer.^2)*dc);
npswf = 12; % minimal number of modes for order 0

for r=noise_levels
  r
  % setting noise strength for predefined relative L2-error level
  log_sigma = log(r)+log(l2_data_norm)-log(sqrt(2))-(lgamma(129.0/2.0)-lgamma(64))-0.5*log(dc);
  sigma = exp(log_sigma);

  % check that noise mapping is correct
  % z = sigma*randn(size(bdwdth_grid));
  % sqrt(dc*trapz(z.^2))/l2_data_norm

  nsamples = 1000;
  z = sigma*randn(nsamples, 128);
  recon_noisy = zeros(nsamples, 256);

  [recon_npswf_noiseless, recon_npswf_grid_noiseless] = reconFromHankelPSWF2D(hankel_data=data_integer, ...
                                                             hankel_grid=bdwdth_grid, ...
                                                             c=c, ...
                                                             order=uint8(order_integer),
                                                             npswf=uint8(npswf), ...
                                                             phi_grid_size=256,
                                                             recon_grid_size=256,
                                                             interp_method = "linear", ...
                                                             rpadd = 1.3,
                                                             verbose_level = 0, ...
                                                             do_plotting = false, ...
                                                             plot2d_size = 64
                                                            );

  for i = 1:nsamples
    i
    data_integer_noisy = data_integer + z(i, :);
    [recon_npswf, recon_npswf_grid] = reconFromHankelPSWF2D(hankel_data=data_integer_noisy, ...
                                                             hankel_grid=bdwdth_grid, ...
                                                             c=c, ...
                                                             order=uint8(order_integer),
                                                             npswf=uint8(npswf), ...
                                                             phi_grid_size=256,
                                                             recon_grid_size=256,
                                                             interp_method = "linear", ...
                                                             rpadd = 1.3,
                                                             verbose_level = 0, ...
                                                             do_plotting = false, ...
                                                             plot2d_size = 64
                                                            );
    recon_noisy(i, :) = recon_npswf(:);
  endfor

  recon_mean = mean(recon_noisy, 1);
  recon_std = std(recon_noisy, 1);

  h=figure('position', [10 10 822 422]);
  image_filename = sprintf("../prolate-2d-3d-experiment/images/errors/errs-bands-ord-%d-modes-%d-noise-%2.2f.png", order_integer, npswf, r*100);

  plot(common_grid, recon_mean, sprintf(";mean, noiselevel %2.1f%%;", 100*r), "linewidth", 1.5);
  ylm = ylim;
  xlm = xlim;
  title(sprintf("95%% error bands for npswf, order=%d, modes=%d", order_integer, npswf), 'position', [mean(xlm) ylm(2)]);
  legend("autoupdate", "off");
  patch([common_grid flip(common_grid)], [recon_mean-2*recon_std flip(recon_mean + 2*recon_std)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
  legend("autoupdate", "on");
  ylim(ylm);
  hold on;
  plot(recon_npswf_grid_noiseless, recon_npswf_noiseless, ";npswf noiseless;", "linewidth", 1.5);
  %hold on;
  plot(recon_npswf_grid, arrayfun(phantom_func, recon_npswf_grid), ";phantom;", "linewidth", 1.5);

  legend("location", "northeastoutside");
  set(gca, 'fontsize', 13);
  saveas(h, image_filename);
endfor

%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Choose number of pswfs via residual minimization in 2D


recon_grid_size = 256;
common_recon_grid = linspace(0., 1., recon_grid_size);
dx = common_recon_grid(2) - common_recon_grid(1);
phantom1d_sampled = arrayfun(phantom_func, common_recon_grid);

r = 0.00; % noise level
order_integer = 0.0;
% DEBUG
phantom_func_sqrt = @(x) phantom_func(x)*sqrt(x);

data_integer = arrayfun(@(t) HankelIntegralTransformF(phantom_func, order_integer, t, 0.0, rsupp), ... % DEBUG
    bdwdth_grid);

% TEST OF 'HankelIntegralTransformD' vs 'HankelIntegralTransformF'
##data_integer_test = arrayfun(@(t) HankelIntegralTransformD(phantom1d_sampled(:), ...
##                                                           common_recon_grid(:), ...
##                                                           order_integer, ...
##                                                           t), bdwdth_grid);
##figure;
##plot(bdwdth_grid, data_integer, ";functional integration;");
##hold on;
##plot(bdwdth_grid, data_integer_test, ";discrete integration;");
##hold off;
##legend("location", "northeastoutside");



dc = bdwdth_grid(2) - bdwdth_grid(1);
l2_data_norm = sqrt(trapz(bdwdth_grid, data_integer.^2));
npswf_max = 20;

% generate additive noise of a given noise level
log_sigma = log(r)+log(l2_data_norm)-log(sqrt(2))-(lgamma(129.0/2.0)-lgamma(64))-0.5*log(dc);
sigma = exp(log_sigma);
z = sigma*randn(1, 128);
data_integer_noisy = data_integer + z(1, :);

l2_residuals = zeros(npswf_max + 1, 1);
cont_recon_npswf_noisy = zeros(21, recon_grid_size);
cont_recon_npswf_noisy_grid = zeros(21, recon_grid_size);
cont_data_recon_reproj = zeros(21, length(bdwdth_grid));

l2_npswf_recon_residuals = zeros(npswf_max + 1, 1);
l2_naive_recon_residual  = 0.;

[recon_naive_noisy, recon_naive_noisy_grid] = reconFromHankelNaive(data_integer_noisy, bdwdth_grid, c, order_integer, 256, do_plotting = true, 128);
data_naive_reproj = arrayfun(@(t) HankelIntegralTransformD(recon_naive_noisy(:), ...
                                                           recon_naive_noisy_grid(:), ...
                                                           order_integer, ...
                                                           t), bdwdth_grid);
l2_naive_residual = sqrt(trapz(bdwdth_grid, (data_naive_reproj - data_integer_noisy).^2));
l2_naive_recon_residual = sqrt(trapz(common_recon_grid, (phantom1d_sampled - recon_naive_noisy).^2));

for npswf = 0:npswf_max
  npswf
  [recon_npswf_noisy, recon_npswf_noisy_grid ] = reconFromHankelPSWF2D(
                                                               data_integer_noisy, ...
                                                               bdwdth_grid, ...
                                                               c, ...
                                                               uint8(order_integer),
                                                               uint8(npswf), ...
                                                               512,
                                                               recon_grid_size,
                                                               "linear", ...
                                                               4.0,
                                                               0, ...
                                                               false, ...
                                                               64
                                                              );
  % save reconstructions
  cont_recon_npswf_noisy(npswf + 1, :) = recon_npswf_noisy;
  cont_recon_npswf_noisy_grid(npswf + 1, :) = recon_npswf_noisy_grid;

  % compute re-projection to Hankel transform
  %plot2dfrom1d(real(recon_npswf_noisy), recon_npswf_noisy_grid, sprintf("order %2.1f, modes=%d", order_integer, npswf));
  data_recon_integer = arrayfun(@(t) HankelIntegralTransformD(recon_npswf_noisy, ...
                                                              recon_npswf_noisy_grid, ...
                                                              order_integer, ...
                                                              t), bdwdth_grid);
  cont_data_recon_reproj(npswf + 1, :) = data_recon_integer(:);

  l2_npswf_recon_residuals(npswf + 1) = sqrt(trapz(common_recon_grid, (phantom1d_sampled - recon_npswf_noisy).^2));
  l2_residuals(npswf + 1) = sqrt(trapz(bdwdth_grid, (data_recon_integer-data_integer_noisy).^2));
endfor
[minval, minidx] = min(l2_residuals, [], 1);


h=figure('position', [10 10 822 422]);
plot(recon_naive_noisy_grid, recon_naive_noisy, ";naive;", "linewidth", 1.5);
hold on;
for idx = 10:10
  plot(cont_recon_npswf_noisy_grid(idx, :), cont_recon_npswf_noisy(idx, :), sprintf(";opt-residual npswf, modes=%d;", idx), "linewidth", 1.5);
endfor
plot(common_recon_grid, arrayfun(phantom_func, common_recon_grid), ";phantom;", "linewidth", 1.5);
legend("location", "northeastoutside");
title(sprintf("order=%2.1f, noiselevel=%2.1f%%", order_integer, 100*r));
set(gca, "fontsize", 14);
hold off;

% plot hankel transforms
h=figure('position', [10 10 822 422]);
plot(bdwdth_grid, data_integer, ";original;", "linewidth", 1.5);
hold on;
plot(bdwdth_grid, data_integer_noisy, ";noisy;", "linewidth", 1.5);
for idx=12:12
  plot(bdwdth_grid, cont_data_recon_reproj(idx, :), sprintf(";opt-residual recon re-projection, modes=%d;", idx), "linewidth", 1.5);
endfor
plot(bdwdth_grid, data_naive_reproj, sprintf(";naive re-projection;"), "linewidth", 1.5);
set(gca, "fontsize", 14);
hold off;


figure;
plot(0:npswf_max, l2_residuals, ";npswf;", "linewidth", 1.5);
hold on;
plot(0:npswf_max, l2_naive_residual*ones(1, npswf_max + 1), ";naive;", "linewidth", 1.5);
ylim([0 1]);
hold off;
title("data l2-residuals");
legend("location", "northeastoutside");
set(gca, "fontsize", 14);

figure;
plot(0:npswf_max, l2_npswf_recon_residuals, ";npswf;");
hold on;
plot(0:npswf_max, l2_naive_recon_residual*ones(1, npswf_max + 1), ";naive;");
ylim([0 1]);
hold off;
title("recon l2-residuals");
legend("location", "northeastoutside");
set(gca, "fontsize", 14);







