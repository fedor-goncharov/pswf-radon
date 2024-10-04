% demonstrate how to scale spatially and by value sinogram before iradon usage

ball_sino_func = @(x) (x >= -1 && x <= 1)*2*sqrt(1-x^2);
ball_sinogram = arrayfun(ball_sino_func, linspace(-1,1, 256));
ball_sinogram = ones(512, 1)*ball_sinogram;


% incorrect padding
ball_recon_test = iradon(ball_sinogram', linspace(0, 180, 512));
figure, imagesc(ball_recon_test, [min(min(ball_recon_test)), max(max(ball_recon_test))]), colorbar;
set(gca, 'fontsize', 13);

% do correct padding
figure, imagesc(ball_sinogram, [min(min(ball_sinogram)), max(max(ball_sinogram))]);
ball_sinogram_padded = padarray(ball_sinogram, [0, 54]); %  53. = (sqrt(2)-1)*128
ball_recon_test = iradon(ball_sinogram_padded', linspace(0, 180, 512));

figure('position', [10 10 400 300]), imagesc(ball_recon_test, [min(min(ball_recon_test)), max(max(ball_recon_test))]);
set(gca, 'fontsize', 13);
colorbar;
ball_recon_test_scaled = ball_recon_test*128;
figure, imagesc(ball_recon_test_scaled, [min(min(ball_recon_test_scaled)), max(max(ball_recon_test_scaled))]), colorbar;
figure, plot(linspace(-1,1, 256), ball_recon_test_scaled(128, :), 'linewidth', 1.5);
xlim([-1.5, 1.5]);
ylim([-0.2, 1.5]);
