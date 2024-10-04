% phantom design by Zajtsev Rodion zajtsev.rv@phystech.edu
% original code is due to Sabinin Grigory gvsabinin@gmail.com
% code has been refactored/rewrittnen by Fedor Goncharov fedor.goncharov.ol@gmail.com


% generate 'ring' phantom
% ------------------------------------------------------------------------------

amplitude = 1.0;
r1 = 0.3;
r2 = 0.5
r3 = -0.50;
r4 = -0.65;
r5 = -0.8
rend = -0.9;
phantom_func = @(x) amplitude*(x > r1 && x <= r2) ...
  + 1.5*(x > r3 && x <= r4) ...
  + amplitude*(x > r5 && x <= rend);

visualize_phantom1d = true;
if visualize_phantom1d
  figure("position", [0 0 512 512]);
  x = linspace(0, 1., 1000);
  y = arrayfun(phantom_func, x);
  plot(x,y)
  title("1d-phantom");
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
  title("phantom-2d");
endif
% ------------------------------------------------------------------------------



% general settings
% ------------------------------------------------------------------------------

rsupp = 1.0; % support in [0, sigma] in the paper
bdwdth_full = 50;   % maximal measured bandwidth
bdwdth_restr = 20; % its restriction
bdwdth_nsize = 128; % 2048
bdwdth_grid = linspace(0, bdwdth_restr, bdwdth_nsize);
c = rsupp*bdwdth_restr;

% ------------------------------------------------------------------------------


% 2D Hankel data generation
% ------------------------------------------------------------------------------
% set order of Hankel transform and generate data
order_integer = 1;
data_integer = arrayfun(@(t) HankelIntegralTransformF(phantom_func, order_integer, t, 0.0, rsupp), ...
  bdwdth_grid);


% check : plot hankel transforms
figure()
plot(bdwdth_grid, data_integer, sprintf(";Hankel transform for %2.2f;", order_integer));
legend("show");
xlabel("t");

% 2D reconstructions
% ------------------------------------------------------------------------------
% naive reconstruction via forward Hankel transform (it works for growing c)
[recon_naive, recon_naive_grid] = reconFromHankelNaive(data_integer, bdwdth_grid, c, order_integer, 256, do_plotting = true, 128);

% reconstruction using PSWFs
[recon_npswf, recon_npswf_grid] = reconFromHankelPSWF2D(hankel_data=data_integer, ...
                                                         hankel_grid=bdwdth_grid, ...
                                                         c=c, ...
                                                         order=uint8(order_integer),
                                                         npswf=uint8(16), ...
                                                         phi_grid_size=384,
                                                         recon_grid_size=256,
                                                         interp_method = "linear", ...
                                                         rpadd = 1.3,
                                                         verbose_level = 1, ...
                                                         do_plotting = true, ...
                                                         plot2d_size = 64
                                                        );
% ------------------------------------------------------------------------------



% 3D Hankel data generation
% ------------------------------------------------------------------------------
verbose_level = 2;
do_plotting = true;

order_hfinteger = 0.5; % set order for half-integer transform
data_hfinteger = arrayfun(@(t) HankelIntegralTransformF(phantom_func, order_hfinteger, t, 0.0, rsupp), ...
  bdwdth_grid);
% ------------------------------------------------------------------------------

% 3D reconstructions
% ------------------------------------------------------------------------------
[recon_naive, recon_naive_grid] = reconFromHankelNaive(data_hfinteger, bdwdth_grid, c, order_hfinteger, 256, do_plotting = true, 128);


n = order_hfinteger - 0.5;
phi_grid_size = 128;
theta_grid_size = 128;
npswf = 14;
recon_grid_size = 128;
[recon, recon_grid] = reconFromHankelPSWF3D(data_hfinteger, ...
                                            bdwdth_grid, ...
                                            c, ...
                                            uint8(n),
                                            uint8(npswf), ...
                                            theta_grid_size,
                                            phi_grid_size,
                                            recon_grid_size,
                                            interp_method = "linear", ...
                                            rpadd = 0.5,
                                            verbose_level = 0, ...
                                            do_plotting = false, ...
                                            plot2d_size = 128
                                           );

figure, plot(recon_grid, recon, ";PSWF reconstruction;");
hold on;
plot(recon_grid, arrayfun(phantom_func, recon_grid), ";phantom;");
hold on;
plot(recon_grid, recon_naive, ";zero extension reconstruction;");
legend("show");
title(sprintf("3D npswf reconstruction c=%2.2f, order=%2.2f, npswf=%d", c, order_hfinteger, npswf));

% ------------------------------------------------------------------------------

