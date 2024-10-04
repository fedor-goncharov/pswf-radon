% ------------------------------------------------------------------------------
% --------------------------- 2D PHANTOM GENERATION ----------------------------

amp1 = 1.0;
ball1_cent_x = 0.19;
ball1_cent_y = 0.25;
ball1_r = 0.13;
ball1_func = @(x,y) amp1*((x-ball1_cent_x)^2 + (y-ball1_cent_y)^2 < ball1_r^2);

amp2 = 1.0;
ball2_cent_x = 0.34;
ball2_cent_y = 0.56;
ball2_r = 0.13;
ball2_func = @(x,y) amp2*((x-ball2_cent_x)^2 + (y-ball2_cent_y)^2 < ball2_r^2);

amp3 = 1.0;
ball3_cent_x = 0.01;
ball3_cent_y = 0.55;
ball3_r = 0.14;
ball3_func = @(x,y) amp3*((x-ball3_cent_x)^2 + (y-ball3_cent_y)^2 < ball3_r^2);

amp_bnd = 1.0;
ball_bnd_r = 1.0;
ball_bnd_func = @(x,y) amp_bnd*(x^2 + y^2 < ball_bnd_r^2);

phantom_2d_func = @(x,y,bnd_flag) (ball1_func(x,y) + ...
                           ball2_func(x,y) + ...
                           ball3_func(x,y)) + ...
                           (bnd_flag == true)*ball_bnd_func(x,y);

visualize_phantom2d = true;
phantom_add_background = false;

phantom_grid_size = 256;
[XX,YY] = meshgrid(linspace(-1,1, phantom_grid_size), linspace(1,-1,phantom_grid_size));
phantom_2d = arrayfun(phantom_2d_func, XX, YY, phantom_add_background);

if visualize_phantom2d
  phantom_2d_image = phantom_2d;
  phantom_2d_image(XX.^2 + YY.^2 > 1.) = -0.5;
  figure("position", [0 0 532 512]);
  imagesc(phantom_2d_image);
  set(gca, "fontsize", 14);
  colormap([1 1 1; jet(256)]);
  caxis([-0.5 2]);
  axis off;
endif


d12 = sqrt((ball1_cent_x-ball2_cent_x)^2 + (ball1_cent_y-ball2_cent_y)^2) ...
    - ball1_r - ball2_r;
d23 = sqrt((ball2_cent_x-ball3_cent_x)^2 + (ball2_cent_y-ball3_cent_y)^2) ...
    - ball2_r - ball3_r;
d13 = sqrt((ball1_cent_x-ball3_cent_x)^2 + (ball1_cent_y-ball3_cent_y)^2) ...
    - ball1_r - ball3_r;

% postsmoothing

% none

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% --------------------- MULTIPOLE EXPANSION ------------------------------------

phantom_grid_x = linspace(-1,1, phantom_grid_size);
phantom_grid_y = linspace(1, -1, phantom_grid_size);
[XXph, YYph] = meshgrid(phantom_grid_x, phantom_grid_y);
rad_grid_size = 256;
phi_grid_size = 256;
%circle_interp_method = "linear";
divide_2pi = true;
grid_pos_o = false;

K = 20; % maximal harmonic index

% container for multipole expansions
multipole_exps_vals = ones(2*K + 1, rad_grid_size);
multipole_exps_grid = ones(1, rad_grid_size);

for order_ind = 1:2*K+1
  order = -K + (order_ind-1);
  printf("multipole expansions : order %d\n", order);
  [rad_proj, rad_proj_grid] = HarmonicProjection2DToRadial1D(phantom_2d, ...
                                                             XXph, ...
                                                             YYph, ...
                                                             rad_grid_size, ...
                                                             order, ...
                                                             phi_grid_size, ...
                                                             "linear", ...
                                                             divide_2pi, ...
                                                             grid_pos_o
                                                             );
  multipole_exps_vals(order_ind, :) = rad_proj;
  if (order_ind == 1)
    multipole_exps_grid = rad_proj_grid;
  endif
endfor


% --------------- check approximation by multipole expansions ------------------


multipole_exp_appr = zeros(phantom_grid_size, phantom_grid_size);
RRph = sqrt(XXph.^2 + YYph.^2);

for order_ind = 1:2*K+1
  order = -K + (order_ind-1);
  printf("multipole expansions recovery : order %d\n", order);
  radial_2d_func = interp1(multipole_exps_grid, multipole_exps_vals(order_ind, :), RRph, 0.0); % extrapolation with 0 outside of radius
  angular_2d_func = (XXph./RRph + 1i*YYph./RRph).^(-order);
  assert(sum(isnan(angular_2d_func(:))) == 0, "NaN from zero division");
  multipole_exp_appr += radial_2d_func .* angular_2d_func;
endfor

visualize_phantom2d_appr = true;
if visualize_phantom2d_appr
  figure("position", [0 0 532 512]);
  multipole_appr_image = real(multipole_exp_appr);
  multipole_appr_image(XX.^2 + YY.^2 > 1.0) = -0.25;
  imagesc(multipole_appr_image);
  %title(sprintf("multipole expansion, Re, K=%d", K));
  set(gca, "fontsize", 14);
  colormap([1 1 1; jet(256)]);
  caxis([-0.25 1.25]);
  axis off;

  figure("position", [0 0 532 512]);
  imagesc(imag(multipole_exp_appr));
  %title(sprintf("multipole expansion, Im, K=%d", K));
  set(gca, "fontsize", 14);
  colormap jet;
  axis off;
endif

% ------------------------------------------------------------------------------


% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% --------------------- HANKEL DATA GENERATION ---------------------------------
% set bandwidth
rsupp = 1.0; % support in [0, sigma] in the paper
bdwdth_restr = 15; % its restriction
bdwdth_nsize = 128; % 2048
bdwdth_grid = linspace(0, bdwdth_restr, bdwdth_nsize);
c = rsupp*bdwdth_restr;
dc = bdwdth_grid(2) - bdwdth_grid(1);

% u_k -> f_k
f_exps_vals = multipole_exps_vals .* sqrt(multipole_exps_grid);

% hankel data container
hankel_data = ones(2*K+1, bdwdth_nsize);

% computation of Hankel data
for order_ind = 1:2*K+1
  order = -K + (order_ind-1);
  hankel_data(order_ind, :) = arrayfun(@(t) HankelIntegralTransformD(f_exps_vals(order_ind, :), ...
                                                                     multipole_exps_grid, ...
                                                                     order, ...
                                                                     t), bdwdth_grid);
endfor

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% ------------------------ ADD NOISE -------------------------------------------


% add gaussian noise
r = 0.00; % global noise level

l2_2_norms_hankel_data = zeros(1, 2*K+1);
l2_2_total_norm_hankel_data = 0.;
for order_ind = 1:2*K+1
  l2_2_norm = trapz(bdwdth_grid, abs(hankel_data(order_ind, :)).^2);
  l2_2_total_norm_hankel_data += l2_2_norm;
  l2_2_norms_hankel_data(order_ind) = l2_2_norm;
endfor
sigma = r*sqrt(l2_2_total_norm_hankel_data)/sqrt((2*K+1)*c)

hankel_data_noisy = hankel_data;
% sample noise for each
for order_ind = 1:2*K+1
  z = sigma*randn(1, bdwdth_nsize);
  hankel_data_noisy(order_ind, :) += z;
endfor

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% ------------------------ SET DATA FOR RECONSTRUCTIONS ------------------------

hankel_data_to_rec = hankel_data;

reconstruct_from_noisy_data = true;
if (reconstruct_from_noisy_data == true)
  hankel_data_to_rec = hankel_data_noisy;
endif

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------


% ------------------------------------------------------------------------------
% ------------------------ NAIVE RECONSTRUCTIONS -------------------------------

% naive reconstructions
f_recon_exps_naive_vals = ones(2*K+1, rad_grid_size);
u_recon_exps_naive_vals = ones(2*K+1, rad_grid_size);

for order_ind = 1:2*K+1
  order = -K + (order_ind-1);
  [f_recon_exps_naive_vals(order_ind, :), f_recon_grid] = reconFromHankelNaive(hankel_data_to_rec(order_ind, :), ...
                                                                               bdwdth_grid, ...
                                                                               c, ...
                                                                               order, ...
                                                                               rad_grid_size, ...
                                                                               false, ...
                                                                               256);
  u_recon_exps_naive_vals(order_ind, :) = f_recon_exps_naive_vals(order_ind, :) ./ sqrt(f_recon_grid);
endfor
u_recon_exps_naive_vals(:, 1) = 0.; % remove nans
assert(sum(isnan(u_recon_exps_naive_vals(:))) == 0, "NaNs in naive reconstructions");

% visual check
recon_naive_appr = zeros(phantom_grid_size, phantom_grid_size);
RRph = sqrt(XXph.^2 + YYph.^2);

for order_ind = 1:2*K+1
  order = -K + (order_ind-1);
  printf("multipole expansions recovery : order %d\n", order);
  radial_2d_func = interp1(f_recon_grid, u_recon_exps_naive_vals(order_ind, :), RRph, 0.0); % extrapolation with 0 outside of radius
  angular_2d_func = (XXph./RRph + 1i*YYph./RRph).^(-order);
  assert(sum(isnan(angular_2d_func(:))) == 0, "NaN from zero division");
  recon_naive_appr += radial_2d_func .* angular_2d_func;
endfor

visualize_recon_naive_appr = true;
if visualize_recon_naive_appr
  figure("position", [0 0 532 512]);
  naive_recon_image = real(recon_naive_appr);
  naive_recon_image(XX.^2 + YY.^2 > 1.0) = -0.25;
  imagesc(naive_recon_image);
  %title(sprintf("naive multipole expansion, Re, K=%d", K));
  set(gca, "fontsize", 14);
  axis off;
  colormap([1 1 1; jet(256)]);
  caxis([-0.25 1.0]);
  %colorbar;

  figure("position", [0 0 600 512]);
  imagesc(imag(recon_naive_appr));
  title(sprintf("naive multipole expansion, Im, K=%d", K));
  set(gca, "fontsize", 14);
  axis off;
  colorbar;
endif

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% ------------------------ PSWF RECONSTRUCTIONS -------------------------------

npswf_max = 20; % maximal number of harmonics

f_recon_exps_pswf_vals = ones(2*K+1, rad_grid_size);
u_recon_exps_pswf_vals = ones(2*K+1, rad_grid_size);


for order_ind = 1:2*K+1
  order = -K + (order_ind-1);
  order
  % containers for residual minimization
  cont_recon_pswf_vals = zeros(npswf_max+1, rad_grid_size);
  cont_recon_pswf_res =  zeros(1, npswf_max+1); % l2-residuals

  for npswf = 0:npswf_max
    npswf
    [recon_pswf_vals, f_recon_grid ] = reconFromHankelPSWF2D(hankel_data_to_rec(order_ind,:), ...
                                                        bdwdth_grid, ...
                                                        c, ...
                                                        uint8(abs(order)), ...
                                                        uint8(npswf), ...
                                                        256,
                                                        rad_grid_size,
                                                        "linear", ...
                                                        2.0, ...
                                                        0, ...
                                                        false, ...
                                                        64);

    if (order < 0)
      recon_pswf_vals = (-1)^abs(order)*recon_pswf_vals;
    endif

    % save reconstructions to local container
    cont_recon_pswf_vals(npswf+1, :) = recon_pswf_vals;
    % compute local l2-residuals
    recon_pswf_data = arrayfun(@(t) HankelIntegralTransformD(recon_pswf_vals, ...
                                                             f_recon_grid, ...
                                                             order, ...
                                                             t), bdwdth_grid);
    % compute residuals
    cont_recon_pswf_res(npswf+1) = sqrt(trapz(bdwdth_grid, abs(recon_pswf_data-hankel_data_to_rec(order_ind,:)).^2));
  endfor

  % save reconstruction with minimal residual
  [minval, minidx] = min(cont_recon_pswf_res);
  f_recon_exps_pswf_vals(order_ind, :) = cont_recon_pswf_vals(minidx, :);
  u_recon_exps_pswf_vals(order_ind, :) = f_recon_exps_pswf_vals(order_ind, :)./sqrt(f_recon_grid);
endfor

% remove NaNs after zero division
u_recon_exps_pswf_vals(isnan(u_recon_exps_pswf_vals)) = 0.;


% visual check

recon_pswf_appr = zeros(phantom_grid_size, phantom_grid_size);
RRph = sqrt(XXph.^2 + YYph.^2);

for order_ind = 1:2*K+1
  order = -K + (order_ind-1);
  printf("multipole expansions recovery : order %d\n", order);
  radial_2d_func = interp1(f_recon_grid, u_recon_exps_pswf_vals(order_ind, :), RRph, 0.0); % extrapolation with 0 outside of radius
  angular_2d_func = (XXph./RRph + 1i*YYph./RRph).^(-order);
  assert(sum(isnan(angular_2d_func(:))) == 0, "NaN from zero division");
  recon_pswf_appr += radial_2d_func .* angular_2d_func;
endfor

visualize_recon_pswf_appr = true;
if visualize_recon_pswf_appr
  recon_pswf_appr_image = real(recon_pswf_appr);
  recon_pswf_appr_image(XX.^2 + YY.^2 > 1.0) = -0.25;
  h = figure("position", [0 0 532 512]);
  imagesc(recon_pswf_appr_image);
  %title(sprintf("PSWF multipole expansion, Re, K=%d", K));
  set(gca, "fontsize", 14);
  colormap([1 1 1; jet(256)]);
  caxis([-0.25 1.0]);
  axis off;
  print(h, "~/Documents/work/pswf-novikov-isaev-2024/images/pswf-recon-k20-0p00.png", "-dpng", "-S532,512");
  %colorbar;

  %figure("position", [0 0 600 512]);
  %imagesc(imag(recon_pswf_appr));
  %title(sprintf("PSWF multipole expansion, Im, K=%d", K));
  %set(gca, "fontsize", 14);
  %axis off;
  %colorbar;
endif


% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------


