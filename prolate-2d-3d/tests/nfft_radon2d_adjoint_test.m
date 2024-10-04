% test 2d radon reconstruction via NUFFT adjoint

% generate a simple sinogram
% characteristic function of a ball of radius 0.5
sino_ball_func = @(x) 2*(abs(x) < 0.5)*sqrt(0.5^2 - x^2);



grid = linspace(-1, 1, 257);
sinogram_slice = arrayfun(sino_ball_func, grid);
nphi = 128;
sinogram = squeeze(repmat(sinogram_slice, [1 1 nphi]));
figure, imshow(sinogram, [min(min(sinogram)), max(max(sinogram))]);
sinogram = sinogram';

start_f = tic();
[nodes, values, jacobian_weights, Fs] = rtft2d(sinogram, 128, 257, 1.0, 2.0, true, 1); % append many zeros
recon2d = nfft_reconstruct_2d(1024, nodes, values, jacobian_weights, Fs, verbose_level = 1);
printf("Elapsed %f seconds\n", toc(start_f));

figure, imshow(real(recon2d), [min(min(recon2d)), max(max(recon2d))]), colorbar;
figure, plot(linspace(-1,1, 1024), recon2d(:, 512));


% phantom shifted from origin

s_grid_size = 129;
nphi = 128;

s_grid = linspace(-1., 1., s_grid_size);
phi_grid = linspace(0, 2*pi, nphi + 1);
phi_grid = phi_grid(1:end-1);
[PP, SS] = ndgrid(phi_grid, s_grid);

r = 0.15;
shift = [0.0 0.4];
sino_ball_shifted_func = @(phi, s) 2*sqrt(max(0, ...
  r^2-(s-dot(shift, [cos(phi) sin(phi)]))^2));

sinogram = arrayfun(sino_ball_shifted_func, PP, SS);

[nodes, values, jacobian_weights, Fs] = rtft2d(sinogram, nphi, s_grid_size, 1.0, 1.0, 1);
recon2d = nfft_reconstruct_2d(192, nodes, values, jacobian_weights, Fs);
figure, imshow(squeeze(recon2d), [0,1]), colorbar, ...
  title(sprintf("shift=[%2.2f %2.2f]", shift(1), shift(2)));
figure, plot(linspace(-1,1, 192), recon2d(:, 117));


% CONCLUSION:
%         1. module 'rtft2d' is (?) calibrated for usage with NUFFT
%         2. nufft_reconstruct_2d must still be calibrated correctly to know ordering of reconstruction
