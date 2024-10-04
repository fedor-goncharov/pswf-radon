% test 2d radon adjoint

% generate a simple analytic singoram and do reconstruction via NFFT
sino_ball_func = @(x) (abs(x) < 0.15)*pi*(0.15^2 - x^2);

s_grid_size = 129;
sino_1d_slice = arrayfun(sino_ball_func, linspace(-1,1, s_grid_size));
nphi = 128; % from (0, 2*pi)
ntheta = 128; % from (0, pi)
sinogram = squeeze(permute(repmat(sino_1d_slice, [1 1 nphi ntheta]), [1 3 4 2]));


[nodes, values, jacobian_weights, Fs] = rtft3d(sinogram, ntheta, nphi, s_grid_size, 1.0, 0.15, true);
recon3d = nfft_reconstruct_3d(192, nodes, values, jacobian_weights, Fs);
recon3d_slice = real(recon3d(:, :, 96));
figure, imshow(recon3d_slice, [min(min(recon3d_slice)), max(max(recon3d_slice))]), colorbar;
figure, plot(linspace(-1,1,192), recon3d_slice(:, 96));

% phantom shifted from origin

s_grid_size = 129;
nphi = 128;
ntheta = 128;

s_grid = linspace(-1., 1., s_grid_size);
phi_grid = linspace(0, 2*pi, nphi + 1)(1:end-1);
u_grid = linspace(-1., 1., ntheta + 2)(2:end-1);
theta_grid = acos(u_grid);
[PP, TT, SS] = meshgrid(phi_grid, theta_grid, s_grid);

r = 0.15;
shift = [0.2 0.4 -0.5];
sino_ball_shifted_func = @(phi, theta, s) pi*max(0, ...
  r^2-(s-dot(shift, [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)]))^2);

sinogram = arrayfun(sino_ball_shifted_func, PP, TT, SS);

[nodes, values, jacobian_weights, Fs] = rtft3d(sinogram, ntheta, nphi, s_grid_size, 1.0, 0.15, true);
recon3d = nfft_reconstruct_3d(192, nodes, values, jacobian_weights, Fs);
recon3d_slice = real(recon3d(48, :, :));
figure, imshow(squeeze(recon3d_slice), [0,1]), colorbar, ...
  title(sprintf("shift=[%2.2f %2.2f %2.2f]", shift(1), shift(2), shift(3)));

% CONCLUSION: image returned by NFFT is ordered as [z, y, x] with monotone indicies that correspond to (-1, 1)






