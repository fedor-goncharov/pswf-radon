% generation of images for article
% 'PSWF-Radon approach to reconstruction from band-limited Hankel transform'

% ------------------------------------------------------------------------------
% --------------------------- PHANTOM GENERATION -------------------------------
amplitude = 1.0;
r1 = 0.15;
r2 = 0.3;
r3 = 0.5;
r4 = 0.75;
r5 = -1.0;
rend = -1.0;
phantom_func = @(x) amplitude*(x > r1 && x <= r2) ...
  + amplitude*(x > r3 && x <= r4) ...
  + amplitude*(x > r5 && x <= rend);

visualize_phantom1d = true;
if visualize_phantom1d
  figure("position", [0 0 532 512]);
  x = linspace(0, 1., 1000);
  y = arrayfun(phantom_func, x);
  plot(x,y, "linewidth", 1.5);
  %title("phantom");
  %xlabel ("x");
  %ylabel ("y=f(x)");
  ylim([-0.25 1.25])
  %filename_plot = sprintf("phantom_1d.png");
  %print(h, [filename_plot], "-dpng", "-S768,512");
endif% set main parameter 'c' for PSWFs

visualize_phantom2d = true;
if visualize_phantom2d
  figure("position", [0 0 532 512]);
  [XX,YY] = meshgrid(linspace(-1,1, 256), linspace(-1,1,256));
  RR = sqrt(XX.^2 + YY.^2);
  phantom_2d = arrayfun(@(r) phantom_func(r), RR);
  imagesc(phantom_2d, [-0.5 1.5]);
  %title("phantom");
  set(gca, "fontsize", 14);
  axis off;
  %filename_plot = sprintf("phantom.png");
  %print(h, [filename_plot], "-dpng", "-S532,512");
endif

% set bandwidth
rsupp = 1.0; % support in [0, sigma] in the paper
bdwdth_restr = 10; % its restriction
bdwdth_nsize = 128; % 2048
bdwdth_grid = linspace(0, bdwdth_restr, bdwdth_nsize);
c = rsupp*bdwdth_restr;
dc = bdwdth_grid(2) - bdwdth_grid(1);
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------


% ----------------------- RESIDUALS FOR NOISE/NOISELESS CASES IN 2D-------------


recon_grid_size = 256;
common_recon_grid = linspace(0., 1., recon_grid_size);
dx = common_recon_grid(2) - common_recon_grid(1);
phantom1d_sampled = arrayfun(phantom_func, common_recon_grid);

order_integer = 0;
npswf_max = 20;

% generate noiseless data
data_integer = arrayfun(@(t) HankelIntegralTransformF(phantom_func, order_integer, t, 0.0, rsupp), ... % DEBUG
    bdwdth_grid);


% add gaussian noise
r = 0.2; % noise level
l2_data_norm = sqrt(trapz(bdwdth_grid, data_integer.^2));
log_sigma = log(r)+log(l2_data_norm)-log(sqrt(2))-(lgamma(129.0/2.0)-lgamma(64))-0.5*log(dc);
sigma = exp(log_sigma);
z = sigma*randn(1, 128);
data_integer_noisy = data_integer + z(1, :);

% set data containers
cont_recon_npswf_noisy = zeros(21, recon_grid_size);
cont_recon_npswf_noisy_grid = zeros(21, recon_grid_size);
cont_data_recon_reproj = zeros(21, length(bdwdth_grid));

l2_rel_error_mode = true;

l2_npswf_data_residuals = zeros(npswf_max + 1, 1);
l2_npswf_recon_residuals = zeros(npswf_max + 1, 1);
l2_naive_recon_residual  = 0.;
l2_naive_data_residual = 0.;

% naive reconstruction from noisy data
[recon_naive_noisy, recon_naive_noisy_grid] = reconFromHankelNaive(data_integer_noisy, bdwdth_grid, c, order_integer, 256, do_plotting = false, 256);
data_naive_reproj = arrayfun(@(t) HankelIntegralTransformD(recon_naive_noisy(:), ...
                                                           recon_naive_noisy_grid(:), ...
                                                           order_integer, ...
                                                           t), bdwdth_grid);
l2_naive_data_residual = sqrt(trapz(bdwdth_grid, abs(data_naive_reproj - data_integer_noisy).^2));
l2_naive_recon_residual = sqrt(trapz(common_recon_grid, abs(phantom1d_sampled - recon_naive_noisy).^2));

if (l2_rel_error_mode == true)
  l2_naive_data_residual /= sqrt(trapz(bdwdth_grid, abs(data_integer_noisy).^2));
  l2_naive_recon_residual /= sqrt(trapz(common_recon_grid, abs(phantom1d_sampled).^2));
endif

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
                                                               2.0,
                                                               0, ...
                                                               false, ...
                                                               64
                                                              );
  % reconstructions to containers
  cont_recon_npswf_noisy(npswf + 1, :) = recon_npswf_noisy;
  cont_recon_npswf_noisy_grid(npswf + 1, :) = recon_npswf_noisy_grid;

  % re-project to Hankel transform
  data_recon_integer = arrayfun(@(t) HankelIntegralTransformD(recon_npswf_noisy, ...
                                                              recon_npswf_noisy_grid, ...
                                                              order_integer, ...
                                                              t), bdwdth_grid);
  cont_data_recon_reproj(npswf + 1, :) = data_recon_integer(:);

  l2_npswf_recon_residuals(npswf + 1) = sqrt(trapz(common_recon_grid, abs(phantom1d_sampled - recon_npswf_noisy).^2));
  l2_npswf_data_residuals(npswf + 1) = sqrt(trapz(bdwdth_grid, abs(data_recon_integer-data_integer_noisy).^2));

  if (l2_rel_error_mode == true)
    l2_npswf_recon_residuals(npswf + 1) /= sqrt(trapz(common_recon_grid, abs(phantom1d_sampled).^2));
    l2_npswf_data_residuals(npswf + 1)  /= sqrt(trapz(bdwdth_grid, abs(data_integer_noisy).^2));
  endif

endfor

% set error mode in filenames
l2_err_mode_abs = "abs";
l2_err_mode_rel = "rel";
l2_err_mode_key = "";

if (l2_rel_error_mode == true)
  l2_err_mode_key = l2_err_mode_rel;
else
  l2_err_mode_key = l2_err_mode_abs;
endif

save_data_flag = true;
if save_data_flag == true
  dir_to_save = "../prolate-2d-3d-experiment/saved-residuals/";

  % pswf data residuals
  filename_data = sprintf("l2_data_residuals_%s_new_pswf_ord%d_noise%d_c%2.0f.txt", ...
                          l2_err_mode_key, ...
                          uint8(order_integer), ...
                          uint8(r*100), ...
                          c);
  fid = fopen([dir_to_save filename_data], "w");
  for i = 1:npswf_max+1
    fprintf(fid, "%d\t%f\n", i-1, l2_npswf_data_residuals(i));
  endfor
  fclose(fid);

  % pswf recon residuals
  filename_data = sprintf("l2_recon_residuals_%s_new_pswf_ord%d_noise%d_c%2.0f.txt", ...
                          l2_err_mode_key, ...
                          uint8(order_integer), ...
                          uint8(r*100), ...
                          c);
  fid = fopen([dir_to_save filename_data], "w");
  for i = 1:npswf_max+1
    fprintf(fid, "%d\t%f\n", i-1, l2_npswf_recon_residuals(i));
  endfor
  fclose(fid);

  % naive data residuals
  filename_data = sprintf("l2_data_residuals_%s_new_naive_ord%d_noise%d_c%2.0f.txt", ...
                          l2_err_mode_key, ...
                          uint8(order_integer), ...
                          uint8(r*100), ...
                          c);
  fid = fopen([dir_to_save filename_data], "w");
  fprintf(fid, "%f\n", l2_naive_data_residual);
  fclose(fid);

  % naive recon residuals
  filename_data = sprintf("l2_recon_residuals_%s_new_naive_ord%d_noise%d_c%2.0f.txt", ...
                          l2_err_mode_key, ...
                          uint8(order_integer), ...
                          uint8(r*100), ...
                          c);
  fid = fopen([dir_to_save filename_data], "w");
  fprintf(fid, "%f\n", l2_naive_recon_residual);
  fclose(fid);
endif

do_plots_flag = true;
close_plots_flag = false;
save_plots_flag = true;

if do_plots_flag == true
  % minimum residual criterion
  [minval, minidx] = min(l2_npswf_data_residuals, [], 1);
  r_plot_lim = 0.95; % plot limit for 1d and 2d plots
  dir_to_save = "../prolate-2d-3d-experiment/saved-images/";

  % 2d plots for optimal PSWFs reconstructions
  plot_grid_size = 256;
  [XX_plot, YY_plot] = meshgrid(linspace(-1., 1., plot_grid_size), linspace(-1., 1., plot_grid_size));
  RR_plot = sqrt(XX_plot.^2 + YY_plot.^2);

  % 2D naive reconstruction
  h = figure("position", [0 0 648 512]);
  recon_naive_2d_plot = interp1(common_recon_grid, recon_naive_noisy, RR_plot, "linear", 0.0);
  recon_naive_2d_plot(RR_plot > r_plot_lim) = 0.0;
  imagesc(recon_naive_2d_plot, [-0.5 1.5]);
  l2_error = l2_naive_recon_residual;
  colorbar("fontsize", 16);
  axis off;
  %title(sprintf("naive reconstruction: l2-error=%1.3f", l2_error));
  set(gca, "fontsize", 16);
  if save_plots_flag == true
    filename_plot = sprintf("2d_new_naive_ord%d_noise%d_c%2.0f.png", uint8(order_integer), uint8(r*100), c);
    print(h, [dir_to_save filename_plot], "-dpng", "-S648,512");
  endif
  if close_plots_flag == true
    close(h);
  endif


  for idx=minidx
    h = figure("position", [0 0 532 512]);
    recon_pswf_2d_plot = interp1(common_recon_grid, real(cont_recon_npswf_noisy(idx, :)), RR_plot, "linear", 0.0);
    recon_pswf_2d_plot(RR_plot > r_plot_lim) = 0.0;
    imagesc(recon_pswf_2d_plot, [-0.5 1.5]);
    axis off;
    %title(sprintf("minimal residual PSWF reconstruction: modes=%d, l2-error=%1.3f", idx-1, l2_error));
    set(gca, "fontsize", 16);
    if save_plots_flag == true
      filename_plot = sprintf("2d_new_pswf_opt%d_ord%d_noise%d_c%2.0f.png", uint8(idx), uint8(order_integer), uint8(r*100), c);
      print(h, [dir_to_save filename_plot], "-dpng", "-S532,512");
    endif
    if close_plots_flag == true
      close(h);
    endif
  endfor

  h=figure('position', [0 0 768 512]);
  plot_recon_grid = common_recon_grid(common_recon_grid < r_plot_lim);
  plot(plot_recon_grid, arrayfun(phantom_func, plot_recon_grid), "color", "k", ";phantom;", "linewidth", 1.5, "linestyle", ":");
  hold on;
  plot(plot_recon_grid, interp1(recon_naive_noisy_grid, recon_naive_noisy, plot_recon_grid, "linear", 0.0), "color", "k", ";naive;", "linewidth", 1.5, "linestyle", "--");
  for idx = minidx
    plot(plot_recon_grid, interp1(cont_recon_npswf_noisy_grid(idx, :), cont_recon_npswf_noisy(idx, :), plot_recon_grid, "linear", 0.0), "color", "k", ...
    sprintf(";PSWF with modes=%d;", idx), "linewidth", 1.5);
  endfor
  legend("off");
  %legend("location", "north", "fontsize", 16);
  %title(sprintf("PSWF reconstruction: c=%2.1f, order=%2.1f, noiselevel=%2.1f%%", c, order_integer, 100*r));
  set(gca, "fontsize", 16);
  ylim([-0.5 1.5]);
  xlim([0.0 r_plot_lim]);
  hold off;
  if save_plots_flag == true
    filename_plot = sprintf("1d_new_recon_mode_dep_ord%d_noise%d_c%2.0f.png", uint8(order_integer), uint8(r*100), c);
    print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
  endif
  if close_plots_flag == true
    close(h);
  endif

  % plot hankel transforms
  h=figure('position', [0 0 768 512]);
  hold on;
  %plot(bdwdth_grid, data_integer, ";original;", "linewidth", 1.5);
  plot(bdwdth_grid, data_integer_noisy, "color", "k", ";noisy;", "linewidth", 1.5, "linestyle", ":"); % data
  hold on;
  for idx=minidx
    plot(bdwdth_grid, cont_data_recon_reproj(idx, :), "color", "k", sprintf(";PSWF with modes=%d;", idx), "linewidth", 1.5, "linestyle", "-");
  endfor
  plot(bdwdth_grid, data_naive_reproj, sprintf(";naive;"), "color", "k", "linewidth", 1.5, "linestyle", "--");
  %legend("location", "north", "fontsize", 16);
  legend("off");
  set(gca, "fontsize", 16);
  hold off;
  if save_plots_flag == true
    filename_plot = sprintf("1d_new_data_mode_dep_ord%d_noise%d_c%2.0f.png", uint8(order_integer), uint8(r*100), c);
    print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
  endif
  if close_plots_flag == true
    close(h);
  endif


  h=figure("position", [0 0 768 512]);
  plot(1:(npswf_max+1), l2_npswf_data_residuals, "color", "k", ";PSWF;", "linewidth", 1.5, "linestyle", "-");
  hold on;
  plot(1:(npswf_max+1), l2_naive_data_residual*ones(1, npswf_max + 1), "color", "k", ";naive;", "linewidth", 1.5, "linestyle", "--");
  ylim([0 1]);
  xlim([1 npswf_max+1]);
  xticks([1 5 10 15 20]);
  hold off;
  legend("off");
  %legend("location", "northeastoutside", "fontsize", 16);
  set(gca, "fontsize", 16);
  if save_plots_flag == true
    filename_plot = sprintf("1d_new_%s_residuals_data_pswf_ord%d_noise%d_c%2.0f.png", ...
                            l2_err_mode_key, ...
                            uint8(order_integer), uint8(r*100), c);
    print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
  endif
  if close_plots_flag == true
    close(h);
  endif


  h=figure("position", [0 0 768 512]);
  plot(1:npswf_max+1, l2_npswf_recon_residuals, "color", "k", ";PSWF;", "linewidth", 1.5, "linestyle", "-");
  hold on;
  plot(1:npswf_max+1, l2_naive_recon_residual*ones(1, npswf_max + 1), "color", "k", ";naive;", "linewidth", 1.5, "linestyle", "--");
  ylim([0 1]);
  xlim([1 npswf_max+1]);
  xticks([1 5 10 15 20]);
  hold off;
  %title("recon l2-residuals");
  %legend("location", "northeastoutside", "fontsize", 16);
  legend("off");
  set(gca, "fontsize", 16);
  if save_plots_flag == true
    filename_plot = sprintf("1d_new_%s_residuals_recon_pswf_ord%d_noise%d_c%2.0f.png", ...
                            l2_err_mode_key, ...
                            uint8(order_integer), uint8(r*100), c);
    print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
  endif
  if close_plots_flag == true
    close(h);
  endif

endif
% ----------------------- RESIDUALS FOR NOISE/NOISELESS CASES IN 3D-------------


recon_grid_size = 256;
common_recon_grid = linspace(0., 1., recon_grid_size);
dx = common_recon_grid(2) - common_recon_grid(1);
phantom1d_sampled = arrayfun(phantom_func, common_recon_grid);

order_hfinteger = 0.5;
npswf_max = 20;

% generate noiseless data
data_hfinteger = arrayfun(@(t) HankelIntegralTransformF(phantom_func, order_hfinteger, t, 0.0, rsupp), ... % DEBUG
    bdwdth_grid);


% add gaussian noise
r = 0.00; % noise level
l2_data_norm = sqrt(trapz(bdwdth_grid, data_hfinteger.^2));
log_sigma = log(r)+log(l2_data_norm)-log(sqrt(2))-(lgamma(129.0/2.0)-lgamma(64))-0.5*log(dc);
sigma = exp(log_sigma);
z = sigma*randn(1, 128);
data_hfinteger_noisy = data_hfinteger + z(1, :);

% set data containers
l2_residuals = zeros(npswf_max + 1, 1);
cont_recon_npswf_noisy = zeros(21, recon_grid_size);
cont_recon_npswf_noisy_grid = zeros(21, recon_grid_size);
cont_data_recon_reproj = zeros(21, length(bdwdth_grid));

% TO RUN BEFORE PLOTTING
l2_rel_error_mode = true;

l2_npswf_recon_residuals = zeros(npswf_max + 1, 1);
l2_naive_recon_residual  = 0.;

% naive reconstruction from noisy data
[recon_naive_noisy, recon_naive_noisy_grid] = reconFromHankelNaive(data_hfinteger_noisy, bdwdth_grid, c, order_hfinteger, 256, do_plotting = true, 256);
data_naive_reproj = arrayfun(@(t) HankelIntegralTransformD(recon_naive_noisy(:), ...
                                                           recon_naive_noisy_grid(:), ...
                                                           order_hfinteger, ...
                                                           t), bdwdth_grid);
l2_naive_residual = sqrt(trapz(bdwdth_grid, (data_naive_reproj - data_hfinteger_noisy).^2));
l2_naive_recon_residual = sqrt(trapz(common_recon_grid, (phantom1d_sampled - recon_naive_noisy).^2));

% TO RUN BEFORE PLOTTING
if (l2_rel_error_mode == true)
  l2_naive_residual /= sqrt(trapz(bdwdth_grid, abs(data_hfinteger_noisy).^2));
  l2_naive_recon_residual /= sqrt(trapz(common_recon_grid, abs(phantom1d_sampled).^2));
endif

for npswf = 0:npswf_max
  npswf
  [recon_npswf_noisy, recon_npswf_noisy_grid ] = reconFromHankelPSWF3D(
                                                               data_hfinteger_noisy, ...
                                                               bdwdth_grid, ...
                                                               c, ...
                                                               uint8(order_hfinteger-0.5),
                                                               uint8(npswf), ...
                                                               256,
                                                               256,
                                                               recon_grid_size,
                                                               "linear", ...
                                                               1.0,
                                                               0, ...
                                                               false, ...
                                                               128
                                                              );
  % reconstructions to containers
  cont_recon_npswf_noisy(npswf + 1, :) = recon_npswf_noisy;
  cont_recon_npswf_noisy_grid(npswf + 1, :) = recon_npswf_noisy_grid;

  % re-project to Hankel transform
  data_recon_integer = arrayfun(@(t) HankelIntegralTransformD(recon_npswf_noisy, ...
                                                              recon_npswf_noisy_grid, ...
                                                              order_hfinteger, ...
                                                              t), bdwdth_grid);
  cont_data_recon_reproj(npswf + 1, :) = data_recon_integer(:);

  l2_npswf_recon_residuals(npswf + 1) = sqrt(trapz(common_recon_grid, (phantom1d_sampled - recon_npswf_noisy).^2));
  l2_residuals(npswf + 1) = sqrt(trapz(bdwdth_grid, (data_recon_integer-data_hfinteger_noisy).^2));

  if (l2_rel_error_mode == true)
    l2_npswf_recon_residuals(npswf + 1) /= sqrt(trapz(common_recon_grid, abs(phantom1d_sampled).^2));
    l2_residuals(npswf + 1)  /= sqrt(trapz(bdwdth_grid, abs(data_integer_noisy).^2));
  endif
endfor

% set error mode in filenames
l2_err_mode_abs = "abs";
l2_err_mode_rel = "rel";
l2_err_mode_key = "";

if (l2_rel_error_mode == true)
  l2_err_mode_key = l2_err_mode_rel;
else
  l2_err_mode_key = l2_err_mode_abs;
endif

save_data_flag = true;
if save_data_flag == true
  dir_to_save = "../prolate-2d-3d-experiment/saved-residuals/";

  % pswf data residuals
  filename_data = sprintf("l2_data_residuals_%s_new_pswf_ord0%d_noise%d_c%2.0f.txt", l2_err_mode_key, ...
                                                                             uint8(order_hfinteger*10), ...
                                                                             uint8(r*100), ...
                                                                             c);
  fid = fopen([dir_to_save filename_data], "w");
  for i = 1:npswf_max+1
    fprintf(fid, "%d\t%f\n", i-1, l2_residuals(i));
  endfor
  fclose(fid);

  % pswf recon residuals
  filename_data = sprintf("l2_recon_residuals_%s_new_pswf_ord0%d_noise%d_c%2.0f.txt", l2_err_mode_key, ...
                                                                             uint8(order_hfinteger*10), ...
                                                                             uint8(r*100), ...
                                                                             c);
  fid = fopen([dir_to_save filename_data], "w");
  for i = 1:npswf_max+1
    fprintf(fid, "%d\t%f\n", i-1, l2_npswf_recon_residuals(i));
  endfor
  fclose(fid);

  % naive data residuals
  filename_data = sprintf("l2_data_residuals_%s_new_naive_ord0%d_noise%d_c%2.0f.txt", l2_err_mode_key, ...
                                                                             uint8(order_hfinteger*10), ...
                                                                             uint8(r*100), ...
                                                                             c);
  fid = fopen([dir_to_save filename_data], "w");
  fprintf(fid, "%f\n", l2_naive_residual);
  fclose(fid);

  % naive recon residuals
  filename_data = sprintf("l2_recon_residuals_%s_new_naive_ord0%d_noise%d_c%2.0f.txt", l2_err_mode_key, ...
                                                                               uint8(order_hfinteger*10), ...
                                                                               uint8(r*100), ...
                                                                               c);
  fid = fopen([dir_to_save filename_data], "w");
  fprintf(fid, "%f\n", l2_naive_recon_residual);
  fclose(fid);
endif

do_plots_flag = true;
close_plots_flag = false;
save_plots_flag = true;

if do_plots_flag == true
  % minimum residual criterion
  [minval, minidx] = min(l2_residuals, [], 1);
  r_plot_lim = 0.95; % plot limit for 1d and 2d plots
  dir_to_save = "../prolate-2d-3d-experiment/saved-images/";

  % 2d plots for optimal PSWFs reconstructions
  plot_grid_size = 256;
  [XX_plot, YY_plot] = meshgrid(linspace(-1., 1., plot_grid_size), linspace(-1., 1., plot_grid_size));
  RR_plot = sqrt(XX_plot.^2 + YY_plot.^2);

  % 2D naive reconstruction
  h = figure("position", [0 0 648 512]);
  recon_naive_2d_plot = interp1(common_recon_grid, recon_naive_noisy, RR_plot, "linear", 0.0);
  recon_naive_2d_plot(RR_plot > r_plot_lim) = 0.0;
  imagesc(recon_naive_2d_plot, [-0.5 1.5]);
  l2_error = l2_naive_recon_residual;
  colorbar("fontsize", 16);
  axis off;
  %title(sprintf("naive reconstruction: l2-error=%1.3f", l2_error));
  set(gca, "fontsize", 16);
  if save_plots_flag == true
    filename_plot = sprintf("2d_new_naive_ord0%d_noise%d_c%2.0f.png", uint8(order_hfinteger*10), uint8(r*100), c);
    print(h, [dir_to_save filename_plot], "-dpng", "-S648,512");
  endif
  if close_plots_flag == true
    close(h);
  endif


  for idx=minidx
    h = figure("position", [0 0 532 512]);
    recon_pswf_2d_plot = interp1(common_recon_grid, real(cont_recon_npswf_noisy(idx, :)), RR_plot, "linear", 0.0);
    recon_pswf_2d_plot(RR_plot > r_plot_lim) = 0.0;
    imagesc(recon_pswf_2d_plot, [-0.5 1.5]);
    axis off;
    %title(sprintf("minimal residual PSWF reconstruction: modes=%d, l2-error=%1.3f", idx-1, l2_error));
    set(gca, "fontsize", 16);
    if save_plots_flag == true
      filename_plot = sprintf("2d_new_pswf_opt%d_ord0%d_noise%d_c%2.0f.png", uint8(idx), uint8(order_hfinteger*10), uint8(r*100), c);
      print(h, [dir_to_save filename_plot], "-dpng", "-S532,512");
    endif
    if close_plots_flag == true
      close(h);
    endif
  endfor

  h=figure('position', [0 0 768 512]);
  plot_recon_grid = common_recon_grid(common_recon_grid < r_plot_lim);
  plot(plot_recon_grid, arrayfun(phantom_func, plot_recon_grid), "color", "k", ";phantom;", "linewidth", 1.5, "linestyle", ":");
  hold on;
  plot(plot_recon_grid, interp1(recon_naive_noisy_grid, recon_naive_noisy, plot_recon_grid, "linear", 0.0), "color", "k", ";naive;", "linewidth", 1.5, "linestyle", "--");
  for idx = minidx
    plot(plot_recon_grid, interp1(cont_recon_npswf_noisy_grid(idx, :), cont_recon_npswf_noisy(idx, :), plot_recon_grid, "linear", 0.0), "color", "k", ...
    sprintf(";PSWF with modes=%d;", idx), "linewidth", 1.5);
  endfor
  legend("off");
  %legend("location", "north", "fontsize", 16);
  %title(sprintf("PSWF reconstruction: c=%2.1f, order=%2.1f, noiselevel=%2.1f%%", c, order_integer, 100*r));
  set(gca, "fontsize", 16);
  ylim([-0.5 1.5]);
  xlim([0.0 r_plot_lim]);
  hold off;
  if save_plots_flag == true
    filename_plot = sprintf("1d_new_recon_mode_dep_ord0%d_noise%d_c%2.0f.png", uint8(order_hfinteger*10), uint8(r*100), c);
    print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
  endif
  if close_plots_flag == true
    close(h);
  endif

  % plot hankel transforms
  h=figure('position', [0 0 768 512]);
  hold on;
  %plot(bdwdth_grid, data_integer, ";original;", "linewidth", 1.5);
  plot(bdwdth_grid, data_hfinteger_noisy, "color", "k", ";noisy;", "linewidth", 1.5, "linestyle", ":"); % data
  hold on;
  for idx=minidx
    plot(bdwdth_grid, cont_data_recon_reproj(idx, :), "color", "k", sprintf(";PSWF with modes=%d;", idx), "linewidth", 1.5, "linestyle", "-");
  endfor
  plot(bdwdth_grid, data_naive_reproj, sprintf(";naive;"), "color", "k", "linewidth", 1.5, "linestyle", "--");
  %legend("location", "north", "fontsize", 16);
  legend("off");
  set(gca, "fontsize", 16);
  hold off;
  if save_plots_flag == true
    filename_plot = sprintf("1d_new_data_mode_dep_ord0%d_noise%d_c%2.0f.png", uint8(order_hfinteger*10), uint8(r*100), c);
    print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
  endif
  if close_plots_flag == true
    close(h);
  endif


  h=figure("position", [0 0 768 512]);
  plot(1:(npswf_max+1), l2_residuals, "color", "k", ";PSWF;", "linewidth", 1.5, "linestyle", "-");
  hold on;
  plot(1:(npswf_max+1), l2_naive_residual*ones(1, npswf_max + 1), "color", "k", ";naive;", "linewidth", 1.5, "linestyle", "--");
  ylim([0 1]);
  xlim([1 npswf_max+1]);
  xticks([1 5 10 15 20]);
  hold off;
  legend("off");
  %legend("location", "northeastoutside", "fontsize", 16);
  set(gca, "fontsize", 16);
  if save_plots_flag == true
    filename_plot = sprintf("1d_new_%s_residuals_data_pswf_ord0%d_noise%d_c%2.0f.png", l2_err_mode_key, uint8(order_hfinteger*10), uint8(r*100), c);
    print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
  endif
  if close_plots_flag == true
    close(h);
  endif


  h=figure("position", [0 0 768 512]);
  plot(1:npswf_max+1, l2_npswf_recon_residuals, "color", "k", ";PSWF;", "linewidth", 1.5, "linestyle", "-");
  hold on;
  plot(1:npswf_max+1, l2_naive_recon_residual*ones(1, npswf_max + 1), "color", "k", ";naive;", "linewidth", 1.5, "linestyle", "--");
  ylim([0 1]);
  xlim([1 npswf_max+1]);
  xticks([1 5 10 15 20]);
  hold off;
  %title("recon l2-residuals");
  %legend("location", "northeastoutside", "fontsize", 16);
  legend("off");
  set(gca, "fontsize", 16);
  if save_plots_flag == true
    filename_plot = sprintf("1d_new_%s_residuals_recon_pswf_ord0%d_noise%d_c%2.0f.png", l2_err_mode_key, uint8(order_hfinteger*10), uint8(r*100), c);
    print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
  endif
  if close_plots_flag == true
    close(h);
  endif

endif
% -----------------------------     2D  ----------------------------------------
% ---------------------------- sine functions  ---------------------------------


% set bandwidth
rsupp = 1.0; % support in [0, sigma] in the paper
bdwdth_restr = 10; % its restriction
bdwdth_nsize = 128; % 2048
bdwdth_grid = linspace(0, bdwdth_restr, bdwdth_nsize);
c = rsupp*bdwdth_restr;
dc = bdwdth_grid(2) - bdwdth_grid(1);

w_grid_size = 20;
w = linspace(c/2, 2*c, w_grid_size); % linspace(c/2, 10*c, w_grid_size);

order_integer = 0.0;
npswf_max = 20;

recon_grid_size = 256;
common_recon_grid = linspace(0., 1., recon_grid_size);
dx = common_recon_grid(2) - common_recon_grid(1);

r = 0.1; % noise level

l2_rel_error_mode = true;

% set data containers
cont_data_noiseless = zeros(w_grid_size, length(bdwdth_grid));
cont_data_noisy     = zeros(w_grid_size, length(bdwdth_grid));

cont_pswf_data_l2_residuals = zeros(w_grid_size, npswf_max + 1);
cont_pswf_recon_l2_residuals = zeros(w_grid_size, npswf_max + 1);

cont_pswf_recon = zeros(w_grid_size, npswf_max + 1, recon_grid_size);
cont_pswf_recon_grid = zeros(w_grid_size, npswf_max + 1, recon_grid_size);
cont_pswf_data_reproj = zeros(w_grid_size, npswf_max + 1, length(bdwdth_grid));

cont_naive_recon = zeros(w_grid_size, recon_grid_size);
cont_naive_recon_grid = zeros(w_grid_size, recon_grid_size);
cont_naive_data_reproj = zeros(w_grid_size, length(bdwdth_grid));

cont_naive_data_l2_residuals = zeros(w_grid_size, 1);
cont_naive_recon_l2_residuals  = zeros(w_grid_size, 1);


for iw = 1:w_grid_size
  iw

  omega = w(iw);
  phantom_func = @(x) sin(omega*x);

  visualize_phantom1d = false;
  if visualize_phantom1d
    figure("position", [0 0 512 512]);
    x = linspace(0, 1., 1000);
    y = arrayfun(phantom_func, x);
    plot(x,y, "linewidth", 1.5);
  endif

  visualize_phantom2d = false;
  if visualize_phantom2d
    figure("position", [0 0 512 512]);
    [XX,YY] = meshgrid(linspace(-1,1, 256), linspace(-1,1,256));
    RR = sqrt(XX.^2 + YY.^2);
    phantom_2d = arrayfun(@(r) phantom_func(r), RR);
    phantom_2d(RR > 1.0) = 0.;
    imagesc(phantom_2d, [-0.5 1.5]);
    set(gca, "fontsize", 14);
    axis off;
  endif
  phantom1d_sampled = arrayfun(phantom_func, common_recon_grid);
% ------------------------------------------------------------------------------

  % generate noiseless data
  data_integer = arrayfun(@(t) HankelIntegralTransformF(phantom_func, order_integer, t, 0.0, rsupp), ... % DEBUG
      bdwdth_grid);
  cont_data_noiseless(iw, :) = data_integer;


  % add gaussian noise
  l2_data_norm = sqrt(trapz(bdwdth_grid, data_integer.^2));
  log_sigma = log(r)+log(l2_data_norm)-log(sqrt(2))-(lgamma(129.0/2.0)-lgamma(64))-0.5*log(dc);
  sigma = exp(log_sigma);
  z = sigma*randn(1, 128);
  data_integer_noisy = data_integer + z(1, :);
  cont_data_noisy(iw, :) = data_integer_noisy;

  % naive reconstruction from noisy data
  [recon_naive_noisy, recon_naive_noisy_grid] = reconFromHankelNaive(data_integer_noisy, bdwdth_grid, c, order_integer, recon_grid_size, do_plotting = false, 256);
  cont_naive_recon(iw, :) = recon_naive_noisy(:);
  cont_naive_recon_grid(iw, :) = recon_naive_noisy_grid;
  cont_naive_recon_l2_residuals(iw) = sqrt(trapz(common_recon_grid, abs(phantom1d_sampled - recon_naive_noisy).^2));

  % reproject naive reconstruction into Hankel space
  data_naive_reproj = arrayfun(@(t) HankelIntegralTransformD(recon_naive_noisy(:), ...
                                                             recon_naive_noisy_grid(:), ...
                                                             order_integer, ...
                                                             t), bdwdth_grid);
  cont_naive_data_reproj(iw, :) = data_naive_reproj;
  cont_naive_data_l2_residuals(iw) = sqrt(trapz(bdwdth_grid, abs(data_naive_reproj - data_integer_noisy).^2));

  if (l2_rel_error_mode == true)
    cont_naive_data_l2_residuals(iw) /= sqrt(trapz(bdwdth_grid, abs(data_integer_noisy).^2));
    cont_naive_recon_l2_residuals(iw) /= sqrt(trapz(common_recon_grid, abs(phantom1d_sampled).^2));
  endif

  for npswf = 0:npswf_max
    npswf
    [recon_pswf_noisy, recon_pswf_noisy_grid ] = reconFromHankelPSWF2D(
                                                                 data_integer_noisy, ...
                                                                 bdwdth_grid, ...
                                                                 c, ...
                                                                 uint8(order_integer),
                                                                 uint8(npswf), ...
                                                                 768,
                                                                 recon_grid_size,
                                                                 "linear", ...
                                                                 3.0,
                                                                 0, ...
                                                                 false, ...
                                                                 64
                                                                );
    % store to containers
    cont_pswf_recon(iw, npswf + 1, :) = recon_pswf_noisy(:);
    cont_pswf_recon_grid(iw, npswf + 1, :) = recon_pswf_noisy_grid(:);

    % re-project PSWF reconstruction to Hankel transform
    data_pswf_reproj = arrayfun(@(t) HankelIntegralTransformD(recon_pswf_noisy, ...
                                                              recon_pswf_noisy_grid, ...
                                                              order_integer, ...
                                                              t), bdwdth_grid);
    cont_pswf_data_reproj(iw, npswf + 1, :) = data_pswf_reproj(:);

    cont_pswf_recon_l2_residuals(iw, npswf + 1) = sqrt(trapz(common_recon_grid, abs(phantom1d_sampled - recon_pswf_noisy).^2));
    cont_pswf_data_l2_residuals(iw, npswf + 1) = sqrt(trapz(bdwdth_grid, abs(data_pswf_reproj-data_integer_noisy).^2));

    if (l2_rel_error_mode == true)
      cont_pswf_recon_l2_residuals(iw, npswf + 1) /= sqrt(trapz(common_recon_grid, abs(phantom1d_sampled).^2));
      cont_pswf_data_l2_residuals(iw, npswf + 1) /= sqrt(trapz(bdwdth_grid, abs(data_integer_noisy).^2));
    endif
  endfor
endfor;

sine_do_plotting = true;
save_plots_flag = true;
close_plots_flag = true;
save_data_flag = true;


% set error mode in filenames
l2_err_mode_abs = "abs";
l2_err_mode_rel = "rel";
l2_err_mode_key = "";

if (l2_rel_error_mode == true)
  l2_err_mode_key = l2_err_mode_rel;
else
  l2_err_mode_key = l2_err_mode_abs;
endif


if save_data_flag == true
  dir_to_save = "../prolate-2d-3d-experiment/saved-residuals/";

  % save id-to-super-resolution map
  filename_data = sprintf("sin_id_to_superres_map_ord%d_noise%d_c%2.0f.txt", uint8(order_integer), uint8(r*100), c);
  fid = fopen([dir_to_save filename_data], "w");
  fprintf(fid, "id\tomega\tsuper-resolution\n");
  for i = 1:w_grid_size
    fprintf(fid, "%d\t%2.2f\t%2.2f\n", i, w(i), w(i)/c);
  endfor
  fclose(fid);


  for data_ind = 1:w_grid_size
    % pswf data residuals
    filename_data = sprintf("l2_sin_id%d_data_%s_residuals_pswf_ord%d_noise%d_c%2.0f.txt", data_ind, l2_err_mode_key, ...
                                                                               uint8(order_integer), ...
                                                                               uint8(r*100), ...
                                                                               c);
    fid = fopen([dir_to_save filename_data], "w");
    for i = 1:npswf_max+1
      fprintf(fid, "%d\t%f\n", i-1, cont_pswf_data_l2_residuals(data_ind, i));
    endfor
    fclose(fid);

    % pswf recon residuals
    filename_data = sprintf("l2_sin_id%d_recon_%s_residuals_pswf_ord%d_noise%d_c%2.0f.txt", data_ind, l2_err_mode_key, ...
                                                                               uint8(order_integer), ...
                                                                               uint8(r*100), ...
                                                                               c);
    fid = fopen([dir_to_save filename_data], "w");
    for i = 1:npswf_max+1
      fprintf(fid, "%d\t%f\n", i-1, cont_pswf_recon_l2_residuals(data_ind, i));
    endfor
    fclose(fid);

    % naive data residuals
    filename_data = sprintf("l2_sin_id%d_data_%s_residuals_naive_ord%d_noise%d_c%2.0f.txt", data_ind, l2_err_mode_key, ...
                                                                               uint8(order_integer), ...
                                                                               uint8(r*100), ...
                                                                               c);
    fid = fopen([dir_to_save filename_data], "w");
    fprintf(fid, "%f\n", cont_naive_data_l2_residuals(data_ind));
    fclose(fid);

    % naive recon residuals
    filename_data = sprintf("l2_sin_id%d_recon_%s_residuals_naive_ord%d_noise%d_c%2.0f.txt", data_ind, l2_err_mode_key, ...
                                                                                 uint8(order_integer), ...
                                                                                 uint8(r*100), ...
                                                                                 c);
    fid = fopen([dir_to_save filename_data], "w");
    fprintf(fid, "%f\n", cont_naive_recon_l2_residuals(data_ind));
    fclose(fid);
  endfor
endif


if (sine_do_plotting == true)
  for plot_ind = 1:w_grid_size
    plot_ind
    dir_to_save = "../prolate-2d-3d-experiment/saved-images/";

    omega = w(plot_ind);
    phantom_func = @(x) sin(omega*x);

    printf("Super-resolution ratio: %2.2f\n", omega/c);

    h=figure("position", [0 0 532 512]);
    [XX,YY] = meshgrid(linspace(-1,1, 256), linspace(-1,1,256));
    RR = sqrt(XX.^2 + YY.^2);
    phantom_2d = arrayfun(@(r) phantom_func(r), RR);
    phantom_2d(RR > 1.0) = 0.;
    imagesc(phantom_2d, [-0.5 1.5]);
    set(gca, "fontsize", 14);
    axis off;
    if save_plots_flag == true
      filename_plot = sprintf("2d_phantom_sin_id%d_ord%d.png", plot_ind, order_integer);
      print(h, [dir_to_save filename_plot], "-dpng", "-S532,512");
    endif
    if close_plots_flag == true
      close(h);
    endif


    % minimum residual criterion
    [minval, minidx] = min(cont_pswf_data_l2_residuals(plot_ind, :), [], 2); % DEBUG (plot_ind, ...)

    % 2d plots for optimal PSWFs reconstructions
    plot_grid_size = 256;
    [XX_plot, YY_plot] = meshgrid(linspace(-1., 1., plot_grid_size), linspace(-1., 1., plot_grid_size));
    RR_plot = sqrt(XX_plot.^2 + YY_plot.^2);

    h=figure("position", [0 0 648 512]);
    recon_naive_2d_plot = interp1(common_recon_grid(:), cont_naive_recon(plot_ind, :)(:), RR_plot, "linear", 0.0);
    imagesc(recon_naive_2d_plot, [-1.5 1.5]);
    colorbar("fontsize", 16);
    axis off;
    %title(sprintf("naive reconstruction: l2-error=%1.3f", l2_error));
    set(gca, "fontsize", 16);
    if save_plots_flag == true
      filename_plot = sprintf("2d_sin_id%d_naive_ord%d_noise%d_c%2.0f.png", plot_ind, uint8(order_integer), uint8(100*r), c);
      print(h, [dir_to_save filename_plot], "-dpng", "-S648,512");
    endif
    if close_plots_flag == true
      close(h);
    endif


    for idx=minidx
      recon_pswf_2d_plot = interp1(common_recon_grid(:), real(cont_pswf_recon(plot_ind, idx, :))(:), RR_plot, "linear", 0.0);
      h=figure("position", [0 0 532 512]);
      imagesc(recon_pswf_2d_plot, [-0.5 1.5]);
      axis off;
      if save_plots_flag == true
          filename_plot = sprintf("2d_sin_id%d_pswf_opt%d_ord%d_noise%d_c%2.0f.png", plot_ind, minidx, uint8(order_integer), uint8(100*r), c);
          print(h, [dir_to_save filename_plot], "-dpng", "-S532,512");
      endif
      if close_plots_flag == true
        close(h);
      endif
    endfor


    % 1d reconstructions
    h=figure('position', [0 0 768 512]);
    plot(common_recon_grid, arrayfun(phantom_func, common_recon_grid), "color", "k", ";phantom;", "linewidth", 1.5, "linestyle", ":");
    hold on;
    plot(common_recon_grid, cont_naive_recon(plot_ind, :), "color", "k", ";naive;", "linewidth", 1.5, "linestyle", "--");
    for idx = minidx
      plot(common_recon_grid, cont_pswf_recon(plot_ind, idx, :), "color", "k", sprintf(";PSWF with modes=%d;", idx), "linewidth", 1.5, "linestyle", "-");
    endfor
    %legend("location", "north", "fontsize", 16);
    legend("off");
    ylim([-1.5 1.5]);
    set(gca, "fontsize", 16);
    hold off;
    if save_plots_flag == true
          filename_plot = sprintf("1d_sin_id%d_recon_mode_dep_ord%d_noise%d_c%2.0f.png", plot_ind, uint8(order_integer), uint8(100*r), c);
          print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
    endif
    if close_plots_flag == true
      close(h);
    endif

    % 1d reprojected hankel transforms
    h=figure('position', [0 0 768 512]);
    plot(bdwdth_grid, cont_data_noisy(plot_ind, :), "color", "k", ";noisy;", "linewidth", 1.5, "linestyle", ":");
    hold on;
    for idx=minidx
      plot(bdwdth_grid, cont_pswf_data_reproj(plot_ind, idx, :), "color", "k", sprintf(";PSWF with modes=%d;", idx), "linewidth", 1.5, "linestyle", "-");
    endfor
    plot(bdwdth_grid, cont_naive_data_reproj(plot_ind, :), "color", "k", sprintf(";naive;"), "linewidth", 1.5, "linestyle", "--");
    %legend("location", "northeast", "fontsize", 16);
    legend("off");
    set(gca, "fontsize", 16);
    hold off;
    if save_plots_flag == true
      filename_plot = sprintf("1d_sin_id%d_data_mode_dep_ord%d_noise%d_c%2.0f.png", plot_ind, uint8(order_integer), uint8(100*r), c);
      print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
    endif
    if close_plots_flag == true
      close(h);
    endif


    h=figure("position", [0 0 768 512]);
    plot(1:npswf_max+1, cont_pswf_data_l2_residuals(plot_ind, :), "color", "k", ";PSWF;", "linewidth", 1.5, "linestyle", "-");
    hold on;
    plot(1:npswf_max+1, "color", "k", cont_naive_data_l2_residuals(plot_ind, :)*ones(1, npswf_max + 1), ";naive;", "linewidth", 1.5, "linestyle", "--");
    ylim([0 1]);
    xlim([1 npswf_max+1]);
    xticks([1 5 10 15 20]);
    hold off;
    %title("data l2-residuals");
    %legend("location", "northeastoutside", "fontsize", 16);
    legend("off");
    set(gca, "fontsize", 16);
    if save_plots_flag == true
      filename_plot = sprintf("1d_sin_id%d_%s_residuals_data_ord%d_noise%d_c%2.0f.png", plot_ind, l2_err_mode_key, uint8(order_integer), uint8(100*r), c);
      print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
    endif
    if close_plots_flag == true
      close(h);
    endif

    h=figure("position", [0 0 768 512]);
    plot(1:npswf_max+1, cont_pswf_recon_l2_residuals(plot_ind, :), "color", "k", ";PSWF;", "linewidth", 1.5, "linestyle", "-");
    hold on;
    plot(1:npswf_max+1, cont_naive_recon_l2_residuals(plot_ind, :)*ones(1, npswf_max + 1), "color", "k", ";naive;", "linewidth", 1.5, "linestyle", "--");
    ylim([0 1]);
    xlim([1 npswf_max+1]);
    xticks([1 5 10 15 20]);
    hold off;
    %title("recon l2-residuals");
    %legend("location", "northeastoutside", "fontsize", 16);
    legend("off");
    set(gca, "fontsize", 16);
    if save_plots_flag == true
      filename_plot = sprintf("1d_sin_id%d_%s_residuals_recon_ord%d_noise%d_c%2.0f.png", plot_ind, l2_err_mode_key, uint8(order_integer), uint8(100*r), c);
      print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
    endif
    if close_plots_flag == true
      close(h);
    endif
  endfor
endif

% ------------------------------------------------------------------------------



% -----------------------------     3D  ----------------------------------------
% ---------------------------- sine functions  ---------------------------------


% set bandwidth
rsupp = 1.0; % support in [0, sigma] in the paper
bdwdth_restr = 10; % its restriction
bdwdth_nsize = 128; % 2048
bdwdth_grid = linspace(0, bdwdth_restr, bdwdth_nsize);
c = rsupp*bdwdth_restr;
dc = bdwdth_grid(2) - bdwdth_grid(1);

w_grid_size = 20;
w = linspace(c/2, 2*c, w_grid_size); % linspace(c/2, 10*c, w_grid_size);

order_hfinteger = 0.5;
npswf_max = 20;

recon_grid_size = 256;
common_recon_grid = linspace(0., 1., recon_grid_size);
dx = common_recon_grid(2) - common_recon_grid(1);

l2_rel_error_mode = true; % l2-relative error or absolute

r = 0.0; % noise level

% set data containers
cont_data_noiseless = zeros(w_grid_sqize, length(bdwdth_grid));
cont_data_noisy     = zeros(w_grid_size, length(bdwdth_grid));

cont_pswf_data_l2_residuals = zeros(w_grid_size, npswf_max + 1);
cont_pswf_recon_l2_residuals = zeros(w_grid_size, npswf_max + 1);

cont_pswf_recon = zeros(w_grid_size, npswf_max + 1, recon_grid_size);
cont_pswf_recon_grid = zeros(w_grid_size, npswf_max + 1, recon_grid_size);
cont_pswf_data_reproj = zeros(w_grid_size, npswf_max + 1, length(bdwdth_grid));

cont_naive_recon = zeros(w_grid_size, recon_grid_size);
cont_naive_recon_grid = zeros(w_grid_size, recon_grid_size);
cont_naive_data_reproj = zeros(w_grid_size, length(bdwdth_grid));

cont_naive_data_l2_residuals = zeros(w_grid_size, 1);
cont_naive_recon_l2_residuals  = zeros(w_grid_size, 1);


for iw = 1:w_grid_size
  iw

  omega = w(iw);
  phantom_func = @(x) sin(omega*x);

  visualize_phantom1d = true;
  if visualize_phantom1d
    figure("position", [0 0 512 512]);
    x = linspace(0, 1., 1000);
    y = arrayfun(phantom_func, x);
    plot(x,y, "linewidth", 1.5);
  endif

  visualize_phantom2d = true;
  if visualize_phantom2d
    figure("position", [0 0 512 512]);
    [XX,YY] = meshgrid(linspace(-1,1, 256), linspace(-1,1,256));
    RR = sqrt(XX.^2 + YY.^2);
    phantom_2d = arrayfun(@(r) phantom_func(r), RR);
    phantom_2d(RR > 1.0) = 0.;
    imagesc(phantom_2d, [-0.5 1.5]);
    set(gca, "fontsize", 14);
    axis off;
  endif
  phantom1d_sampled = arrayfun(phantom_func, common_recon_grid);
% ------------------------------------------------------------------------------

  % generate noiseless data
  data_hfinteger = arrayfun(@(t) HankelIntegralTransformF(phantom_func, order_hfinteger, t, 0.0, rsupp), ... % DEBUG
      bdwdth_grid);
  cont_data_noiseless(iw, :) = data_hfinteger;


  % add gaussian noise
  l2_data_norm = sqrt(trapz(bdwdth_grid, data_hfinteger.^2));
  log_sigma = log(r)+log(l2_data_norm)-log(sqrt(2))-(lgamma(129.0/2.0)-lgamma(64))-0.5*log(dc);
  sigma = exp(log_sigma);
  z = sigma*randn(1, 128);
  data_hfinteger_noisy = data_hfinteger + z(1, :);
  cont_data_noisy(iw, :) = data_hfinteger_noisy;

  % naive reconstruction from noisy data
  [recon_naive_noisy, recon_naive_noisy_grid] = reconFromHankelNaive(data_hfinteger_noisy, bdwdth_grid, c, order_hfinteger, recon_grid_size, do_plotting = true, 256);
  cont_naive_recon(iw, :) = recon_naive_noisy(:);
  cont_naive_recon_grid(iw, :) = recon_naive_noisy_grid;
  cont_naive_recon_l2_residuals(iw) = sqrt(trapz(common_recon_grid, abs(phantom1d_sampled - recon_naive_noisy).^2));

  % reproject naive reconstruction into Hankel space
  data_naive_reproj = arrayfun(@(t) HankelIntegralTransformD(recon_naive_noisy(:), ...
                                                             recon_naive_noisy_grid(:), ...
                                                             order_hfinteger, ...
                                                             t), bdwdth_grid);
  cont_naive_data_reproj(iw, :) = data_naive_reproj;
  cont_naive_data_l2_residuals(iw) = sqrt(trapz(bdwdth_grid, abs(data_naive_reproj - data_hfinteger_noisy).^2));

  for npswf = 0:npswf_max
    npswf
    [recon_pswf_noisy, recon_pswf_noisy_grid ] = reconFromHankelPSWF2D(
                                                                 data_hfinteger_noisy, ...
                                                                 bdwdth_grid, ...
                                                                 c, ...
                                                                 uint8(order_integer),
                                                                 uint8(npswf), ...
                                                                 768,
                                                                 recon_grid_size,
                                                                 "linear", ...
                                                                 3.0,
                                                                 0, ...
                                                                 false, ...
                                                                 64
                                                                );
    % store to containers
    cont_pswf_recon(iw, npswf + 1, :) = recon_pswf_noisy(:);
    cont_pswf_recon_grid(iw, npswf + 1, :) = recon_pswf_noisy_grid(:);

    % re-project PSWF reconstruction to Hankel transform
    data_pswf_reproj = arrayfun(@(t) HankelIntegralTransformD(recon_pswf_noisy, ...
                                                              recon_pswf_noisy_grid, ...
                                                              order_hfinteger, ...
                                                              t), bdwdth_grid);
    cont_pswf_data_reproj(iw, npswf + 1, :) = data_pswf_reproj(:);

    cont_pswf_recon_l2_residuals(iw, npswf + 1) = sqrt(trapz(common_recon_grid, abs(phantom1d_sampled - recon_pswf_noisy).^2));
    cont_pswf_data_l2_residuals(iw, npswf + 1) = sqrt(trapz(bdwdth_grid, abs(data_pswf_reproj-data_hfinteger_noisy).^2));
  endfor
endfor;

sine_do_plotting = true;
save_plots_flag = true;
close_plots_flag = true;
save_data_flag = true;

if save_data_flag == true
  dir_to_save = "../prolate-2d-3d-experiment/saved-residuals/";

  % save id-to-super-resolution map
  filename_data = sprintf("sin_id_to_superres_map_ord0%d_noise%d_c%2.0f.txt", uint8(order_hfinteger*10), uint8(r*100), c);
  fid = fopen([dir_to_save filename_data], "w");
  fprintf(fid, "id\tomega\tsuper-resolution\n");
  for i = 1:w_grid_size
    fprintf(fid, "%d\t%2.2f\t%2.2f\n", i, w(i), w(i)/c);
  endfor
  fclose(fid);


  for data_ind = 1:w_grid_size
    % pswf data residuals
    filename_data = sprintf("l2_sin_id%d_data_residuals_pswf_ord0%d_noise%d_c%2.0f.txt", data_ind,
                                                                               uint8(order_hfinteger*10), ...
                                                                               uint8(r*100), ...
                                                                               c);
    fid = fopen([dir_to_save filename_data], "w");
    for i = 1:npswf_max+1
      fprintf(fid, "%d\t%f\n", i-1, cont_pswf_data_l2_residuals(data_ind, i));
    endfor
    fclose(fid);

    % pswf recon residuals
    filename_data = sprintf("l2_sin_id%d_recon_residuals_pswf_ord0%d_noise%d_c%2.0f.txt", data_ind, ...
                                                                               uint8(order_hfinteger*10), ...
                                                                               uint8(r*100), ...
                                                                               c);
    fid = fopen([dir_to_save filename_data], "w");
    for i = 1:npswf_max+1
      fprintf(fid, "%d\t%f\n", i-1, cont_pswf_recon_l2_residuals(data_ind, i));
    endfor
    fclose(fid);

    % naive data residuals
    filename_data = sprintf("l2_sin_id%d_data_residuals_naive_ord0%d_noise%d_c%2.0f.txt", data_ind, ...
                                                                               uint8(order_hfinteger*10), ...
                                                                               uint8(r*100), ...
                                                                               c);
    fid = fopen([dir_to_save filename_data], "w");
    fprintf(fid, "%f\n", cont_naive_data_l2_residuals(data_ind));
    fclose(fid);

    % naive recon residuals
    filename_data = sprintf("l2_sin_id%d_recon_residuals_naive_ord0%d_noise%d_c%2.0f.txt", data_ind, ...
                                                                                 uint8(order_hfinteger*10), ...
                                                                                 uint8(r*100), ...
                                                                                 c);
    fid = fopen([dir_to_save filename_data], "w");
    fprintf(fid, "%f\n", cont_naive_recon_l2_residuals(data_ind));
    fclose(fid);
  endfor
endif


if (sine_do_plotting == true)

  for plot_ind = 1:w_grid_size
    %plot_ind = 10;

    plot_ind
    dir_to_save = "../prolate-2d-3d-experiment/saved-images/";

    omega = w(plot_ind);
    phantom_func = @(x) sin(omega*x);

    printf("Super-resolution ratio: %2.2f\n", omega/c);

    h=figure("position", [0 0 532 512]);
    [XX,YY] = meshgrid(linspace(-1,1, 256), linspace(-1,1,256));
    RR = sqrt(XX.^2 + YY.^2);
    phantom_2d = arrayfun(@(r) phantom_func(r), RR);
    phantom_2d(RR > 1.0) = 0.;
    imagesc(phantom_2d, [-0.5 1.5]);
    set(gca, "fontsize", 14);
    axis off;
    if save_plots_flag == true
      filename_plot = sprintf("2d_phantom_sin_id%d_ord0%d_noise%d_c%2.0f.png", plot_ind, uint8(order_hfinteger*10), uint8(r*100), c);
      print(h, [dir_to_save filename_plot], "-dpng", "-S532,512");
    endif
    if close_plots_flag == true
      close(h);
    endif


    % minimum residual criterion
    [minval, minidx] = min(cont_pswf_data_l2_residuals(plot_ind, :), [], 2); % DEBUG (plot_ind, ...)

    % 2d plots for optimal PSWFs reconstructions
    plot_grid_size = 256;
    [XX_plot, YY_plot] = meshgrid(linspace(-1., 1., plot_grid_size), linspace(-1., 1., plot_grid_size));
    RR_plot = sqrt(XX_plot.^2 + YY_plot.^2);

    h=figure("position", [0 0 648 512]);
    recon_naive_2d_plot = interp1(common_recon_grid(:), cont_naive_recon(plot_ind, :)(:), RR_plot, "linear", 0.0);
    imagesc(recon_naive_2d_plot, [-1.5 1.5]);
    colorbar("fontsize", 16);
    axis off;
    %title(sprintf("naive reconstruction: l2-error=%1.3f", l2_error));
    set(gca, "fontsize", 16);
    if save_plots_flag == true
      filename_plot = sprintf("2d_sin_id%d_naive_ord0%d_noise%d_c%2.0f.png", plot_ind, uint8(order_hfinteger*10), uint8(100*r), c);
      print(h, [dir_to_save filename_plot], "-dpng", "-S648,512");
    endif
    if close_plots_flag == true
      close(h);
    endif


    for idx=minidx
      recon_pswf_2d_plot = interp1(common_recon_grid(:), real(cont_pswf_recon(plot_ind, idx, :))(:), RR_plot, "linear", 0.0);
      h=figure("position", [0 0 532 512]);
      imagesc(recon_pswf_2d_plot, [-1.5 1.5]);
      axis off;
      if save_plots_flag == true
          filename_plot = sprintf("2d_sin_id%d_pswf_opt%d_ord0%d_noise%d_c%2.0f.png", plot_ind, minidx, uint8(order_hfinteger), uint8(100*r), c);
          print(h, [dir_to_save filename_plot], "-dpng", "-S532,512");
      endif
      if close_plots_flag == true
        close(h);
      endif
    endfor


    % 1d reconstructions
    h=figure('position', [0 0 768 512]);
    plot(common_recon_grid, arrayfun(phantom_func, common_recon_grid), "color", "k", ";phantom;", "linewidth", 1.5, "linestyle", ":");
    hold on;
    plot(common_recon_grid, cont_naive_recon(plot_ind, :), "color", "k", ";naive;", "linewidth", 1.5, "linestyle", "--");
    for idx = minidx
      plot(common_recon_grid, cont_pswf_recon(plot_ind, idx, :), "color", "k", sprintf(";PSWF with modes=%d;", idx), "linewidth", 1.5, "linestyle", "-");
    endfor
    %legend("location", "north", "fontsize", 16);
    legend("off");
    ylim([-1.5 1.5]);
    set(gca, "fontsize", 16);
    hold off;
    if save_plots_flag == true
          filename_plot = sprintf("1d_sin_id%d_recon_mode_dep_ord0%d_noise%d_c%2.0f.png", plot_ind, uint8(order_hfinteger*10), uint8(100*r), c);
          print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
    endif
    if close_plots_flag == true
      close(h);
    endif

    % 1d reprojected hankel transforms
    h=figure('position', [0 0 768 512]);
    plot(bdwdth_grid, cont_data_noisy(plot_ind, :), "color", "k", ";noisy;", "linewidth", 1.5, "linestyle", ":");
    hold on;
    for idx=minidx
      plot(bdwdth_grid, cont_pswf_data_reproj(plot_ind, idx, :), "color", "k", sprintf(";PSWF with modes=%d;", idx), "linewidth", 1.5, "linestyle", "-");
    endfor
    plot(bdwdth_grid, cont_naive_data_reproj(plot_ind, :), "color", "k", sprintf(";naive;"), "linewidth", 1.5, "linestyle", "--");
    %legend("location", "northeast", "fontsize", 16);
    legend("off");
    set(gca, "fontsize", 16);
    hold off;
    if save_plots_flag == true
      filename_plot = sprintf("1d_sin_id%d_data_mode_dep_ord0%d_noise%d_c%2.0f.png", plot_ind, uint8(order_hfinteger*10), uint8(100*r), c);
      print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
    endif
    if close_plots_flag == true
      close(h);
    endif


    h=figure("position", [0 0 768 512]);
    plot(1:npswf_max+1, cont_pswf_data_l2_residuals(plot_ind, :), "color", "k", ";PSWF;", "linewidth", 1.5, "linestyle", "-");
    hold on;
    plot(1:npswf_max+1, "color", "k", cont_naive_data_l2_residuals(plot_ind, :)*ones(1, npswf_max + 1), ";naive;", "linewidth", 1.5, "linestyle", "--");
    ylim([0 1]);
    xlim([1 npswf_max+1]);
    xticks([1 5 10 15 20]);
    hold off;
    %title("data l2-residuals");
    %legend("location", "northeastoutside", "fontsize", 16);
    legend("off");
    set(gca, "fontsize", 16);
    if save_plots_flag == true
      filename_plot = sprintf("1d_sin_id%d_residuals_data_ord0%d_noise%d_c%2.0f.png", plot_ind, uint8(order_hfinteger*10), uint8(100*r), c);
      print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
    endif
    if close_plots_flag == true
      close(h);
    endif

    h=figure("position", [0 0 768 512]);
    plot(1:npswf_max+1, cont_pswf_recon_l2_residuals(plot_ind, :), "color", "k", ";PSWF;", "linewidth", 1.5, "linestyle", "-");
    hold on;
    plot(1:npswf_max+1, cont_naive_recon_l2_residuals(plot_ind, :)*ones(1, npswf_max + 1), "color", "k", ";naive;", "linewidth", 1.5, "linestyle", "--");
    ylim([0 1]);
    xlim([1 npswf_max+1]);
    xticks([1 5 10 15 20]);
    hold off;
    %title("recon l2-residuals");
    %legend("location", "northeastoutside", "fontsize", 16);
    legend("off");
    set(gca, "fontsize", 16);
    if save_plots_flag == true
      filename_plot = sprintf("1d_sin_id%d_residuals_recon_ord0%d_noise%d_c%2.0f.png", plot_ind, uint8(order_hfinteger*10), uint8(100*r), c);
      print(h, [dir_to_save filename_plot], "-dpng", "-S768,512");
    endif
    if close_plots_flag == true
      close(h);
    endif
  endfor
endif

% ------------------------------------------------------------------------------



