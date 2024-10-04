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

% construct 2d-phantom
visualize_phantom2d = true;
phantom2d_size = 128;
linx = linspace(-1.,1., phantom2d_size + 1)(1:end-1);
liny = linspace(-1.,1., phantom2d_size + 1)(1:end-1);
[XX,YY] = meshgrid(linx, liny);
RR = sqrt(XX.^2 + YY.^2);
phantom_2d = arrayfun(@(r) phantom_func(r), RR);

if visualize_phantom2d
  figure;
  imshow(phantom_2d, [min(min(phantom_2d)), max(max(phantom_2d))]);
  title("phantom-2d in experiment with Gaussian noise");
endif

% scale phantom for Fourier space
phantom_2d_farg = phantom_2d ./ sqrt(RR);
phantom_2d_farg(isnan(phantom_2d_farg)) = 0.0;
phantom_2d_farg_padd = padarray(phantom_2d_farg, 4*[phantom2d_size phantom2d_size]);
figure, imagesc(phantom_2d_farg_padd); % plot phantom


% take its Fourier transform
dx = linx(2) - linx(1);
dy = dx;
phantom_2d_f = fftshift(fft2(ifftshift(phantom_2d_farg_padd)))*dx*dy;
figure, imagesc(real(phantom_2d_f)), colorbar;


% Fourier grid for Fourier transform of phantom
Fs = 1./dx;
freqs_lin = linspace(-Fs/2, Fs/2, size(phantom_2d_f, 1) + 1)(1:end-1);
[freqs_XX, freqs_YY] = meshgrid(freqs_lin, freqs_lin);
freqs_sqrt_RR = sqrt(sqrt(freqs_XX.^2 + freqs_YY.^2));


% spherical averaging of Fourier data (should look like Hankel transform of 1D-phantom)
phi_grid_size = 256;
phi_grid = linspace(0, 2*pi, phi_grid_size + 1)(1:end-1);
dphi = phi_grid(2) - phi_grid(1);

xs = vec(cos(phi_grid));
ys = vec(sin(phi_grid));

fourier_r_grid_size = 256;
c=10;
fourier_r_grid = linspace(0, c/(2*pi), fourier_r_grid_size);

% generate samples
sigma = 0.0
nsamples = 1;
hankel_container = zeros(nsamples, fourier_r_grid_size);
hankel_r_grid = fourier_r_grid*2*pi; % re-grid frequencies to canonical hankel transform



for n = 1:nsamples
  n
  noise_f = sigma*randn(size(phantom_2d_f));
  phantom_2d_fn = phantom_2d_f + noise_f;
  phantom_2d_fns = phantom_2d_fn .* freqs_sqrt_RR;
  fourier_sph_means = arrayfun(@(r) trapz(interp2(freqs_XX, freqs_YY, phantom_2d_fns, r*xs, r*ys))*dphi, fourier_r_grid);

  hankel_data_from_fourier = fourier_sph_means / (2*pi)^(1.5);
  hankel_container(n, :) = hankel_data_from_fourier(:);
endfor

% plot mean
figure, plot(hankel_r_grid, mean(real(hankel_container), 1));
title("mean H_0");
set(gca, "fontsize", 14);

% plot standard deviation (note that for 1 sample std is not defined)
figure, plot(hankel_r_grid, std(real(hankel_container), 1));
title("std H_0");
set(gca, "fontsize", 14);


% plot autocovariance
figure, imagesc(cov(real(hankel_container)));
title("autocovariance H_0");
set(gca, "fontsize", 14);


% CHECK AGAINST CANONICAL HANKEL DATA ------------------------------------------

% canonical hankel transform vs. from fourier
hankel_order = 0;
rsupp = 1.0;
bdwdth_grid_size = 128;
bdwdth_grid = linspace(0, c, bdwdth_grid_size);
hankel_data_canonical = arrayfun(@(t) HankelIntegralTransformF(phantom_func, hankel_order, t, 0.0, rsupp), bdwdth_grid);

%dc = bdwdth_grid(2)-bdwdth_grid(1);
%sqrt(trapz(hankel_data.^2)*dc)

figure, plot(bdwdth_grid, hankel_data_canonical, ";canonical;");
hold on;
plot(bdwdth_grid, real(interp1(hankel_r_grid, real(hankel_container(1, :)), bdwdth_grid, "linear", 0.0)), ";from fourier;");
hold off;
legend("location", "northeastoutside");


% NOISE TEST -------------------------------------------------------------------

% reconstruct signal on [0,1] and compute residual in 2D Fourier space
% nothing should change: l2 distance in Fourier space over the ball equals L2-distance
%                        between Hankel transforms

hankel_order_integer = 0;
hankel_data = real(hankel_container(1, :)); % cast to real
hankel_grid = hankel_r_grid;
recon_grid_size = 256;

% set containers
l2_fourier_2d_residuals = zeros(21, 1);
hankel_recons = zeros(21, recon_grid_size);

% naive
[recon_naive, recon_naive_grid] = reconFromHankelNaive(hankel_data, hankel_grid, c, order_integer, 256, do_plotting = false, 128);
hankel_naive_reproj = arrayfun(@(t) HankelIntegralTransformD(recon_naive(:), ...
                                                           recon_naive_grid(:), ...
                                                           hankel_order_integer, ...
                                                           t), hankel_grid);

hankel_naive_residual = sqrt(trapz(hankel_grid, (hankel_naive_reproj(:) - hankel_data(:)).^2));



for npswf = 0:20

  npswf
  % DEBUG
  % reconstruct using pswfs
  [recon_npswf, recon_npswf_grid] = reconFromHankelPSWF2D(hankel_data, ...
                                                          hankel_grid, ...
                                                          c, ...
                                                          uint8(hankel_order_integer),
                                                          uint8(npswf), ...
                                                          256, ...
                                                          recon_grid_size, ...
                                                          "linear", ...
                                                          1.0, ...
                                                          0, ...
                                                          false, ...
                                                          64
                                                          );
  hankel_recons(npswf + 1, :) = real(recon_npswf);

  % cast to real
  recon_npswf = real(recon_npswf);

  % cast reconstruction to 2D Fourier data
  recon_npswf_2d = interp1(recon_npswf_grid, recon_npswf, RR, "linear", 0.0); % IMPORTANT : RR is taken from phantom - grid is one of phantom now on!

  % padd with zeros as original phantom and compute Fourier transform
  recon_npswf_2d(RR > 1.0) = 0.0;
  recon_npswf_2d_farg = recon_npswf_2d ./ sqrt(RR); % DEBUG
  recon_npswf_2d_farg(isnan(recon_npswf_2d_farg)) = 0.0;
  recon_npswf_2d_farg_padd = padarray(recon_npswf_2d_farg, 4*[phantom2d_size phantom2d_size]);
  recon_npswf_2d_f = fftshift(fft2(ifftshift(recon_npswf_2d_farg_padd)))*dx*dy; % fourier transform

  freqs_RR = sqrt(freqs_XX.^2 + freqs_YY.^2);
  freqs_mask_c = zeros(size(freqs_RR));
  freqs_mask_c(freqs_RR < c/(2*pi)) = 1.0;

  % l2-norm between fourier balls of noisy data and the re-projected Fourier data
  df = freqs_lin(2) - freqs_lin(1);
  l2_fourier_2d_residual = sqrt(sum((real((recon_npswf_2d_f - phantom_2d_f).*freqs_mask_c).^2)(:))*df*df/(2*pi));
  l2_fourier_2d_residuals(npswf + 1) = l2_fourier_2d_residual;

endfor

% compute canonical residuals
hankel_canonical_residuals = zeros(21, 1);
dc = hankel_grid(2)-hankel_grid(1);

plot(hankel_grid, hankel_data, ";original data;", "linewidth", 1.5);
hold on;
plot(hankel_grid, hankel_naive_reproj, ";naive;", "linewidth", 1.5);

for npswf = 0:1:15
  hankel_canonical_reproj = arrayfun(@(t) HankelIntegralTransformD(
                                                              hankel_recons(npswf + 1, :), ...
                                                              recon_npswf_grid, ...
                                                              hankel_order, ...
                                                              t), hankel_r_grid);

  plot(hankel_r_grid, hankel_canonical_reproj, sprintf(";pswf modes=%d;", npswf), "linewidth", 1.5);
  hankel_canonical_residuals(npswf + 1) = sqrt(trapz(hankel_r_grid, (hankel_canonical_reproj - hankel_data).^2));

endfor
title("Hankel data through PSWFs re-projections; order 0;");
hold off;

% DON'T FORGET TO RUN THE PLOTTING ABOVE
figure;
plot(0:13, ones(14,1)*hankel_naive_residual, ";naive;", "linewidth", 1.5);
hold on;
plot(0:13, l2_fourier_2d_residuals(1:14), ";through fourier;");
hold on;
plot(0:13, hankel_canonical_residuals(1:14), ";through HankelD;");
hold off;
title("l2-residuals in 2D Fourier space");
legend("location", "northeastoutside");




% -----check Hankel transform L2-norm vs L2-Fourier norms correspondances ------

##dc = hankel_grid(2)-hankel_grid(1);
##hankel_data_l2 = sqrt(trapz(hankel_data.^2)*dc);
##
##freqs_RR = sqrt(freqs_XX.^2 + freqs_YY.^2);
##freqs_mask_c = zeros(size(freqs_RR));
##freqs_mask_c(freqs_RR < c/(2*pi)) = 1.0;
##df = freqs_lin(2)-freqs_lin(1);
##
##% integration formula is a bit poor as well the mask
##hankel_data_l2_from_fourier = sqrt(sum((real(phantom_2d_f.*freqs_mask_c).^2)(:))*df*df/(2*pi));
##
##hankel_data_l2
##hankel_data_l2_from_fourier

