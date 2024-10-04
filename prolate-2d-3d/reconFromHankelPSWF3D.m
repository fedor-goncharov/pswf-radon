% reconstruction from Hankel transform of half-integer order (order = n + 1/2) limited on [0, c]
% via PSWFs SVD decomposition and Radon inversion in 3D

function [recon, recon_grid] = reconFromHankelPSWF3D(hankel_data, ...
                                                 hankel_grid, ...
                                                 c, ...
                                                 n,
                                                 npswf, ...
                                                 theta_grid_size,
                                                 phi_grid_size,
                                                 recon_grid_size,
                                                 interp_method = "linear", ...
                                                 rpadd = 0.5,
                                                 verbose_level = 0, ...
                                                 do_plotting = false, ...
                                                 plot2d_size = 64
                                                 )
  % Description:
  %   ...
  % params:
  %   hankel_data         :
  %   hankel_grid         :
  %   c                   :
  %   n                   :
  %   npswf               :
  %   phi_grid_size       :
  %   theta_grid_size     :
  %   recon_grid_size     :
  %   interp_method       :
  %   rpadd               :
  %   verbose_level       :
  %   do_plotting         :
  %   plot2d_size         :
  %
  % assumed:
  %   reconstructed signal is supported on [0, 1]

  % input checks
  assert(hankel_grid(1) == 0.0 && hankel_grid(end) == c, "Problem with grid - check doc for correct usage. Abort. \n");
  assert(isinteger(npswf) && npswf >= uint8(0), "Provided 'npswf' is not of integer type or negative (e.g. npswf=uint8(X)). Abort.\n");
  assert(isinteger(n), "Provided order is not of integer type. Did you cast it to uintX(Y)?\n");
  assert(n >= uint8(0), "Provided order is negative. Abort.\n")

  n = double(n);
  order = n + 0.5;

  % global timer
  start = tic;

  if (verbose_level >= 1)
    printf("Reconstruction from Hankel transform data of signal supported in [0,1]\n");
    printf("Main parameters\n");
    printf("\tData length: %d, grid start: %f, grid end: %f\n", length(hankel_data), hankel_grid(1), hankel_grid(end));
    printf("\tBandwidth 'c': %f\n", c);
    printf("\tHankel transform order: %f\n", order);
    printf("\tNumber of PSWFs used: %d\n", npswf);
    printf("\tNumber of projections for Radon embedding: %d\n", phi_grid_size*theta_grid_size);
    printf("\tPad Radon data with zeros up to: %f\n", rpadd);
    printf("Auxiliary parameters\n");
    printf("\tInterpolation method for angle averaging: %s\n", interp_method);
    printf("\tDo plots: %d\n", do_plotting);
  endif

  % plot data before
  if (do_plotting)
    figure, ...
    plot(hankel_grid, hankel_data), ...
    title(sprintf("Hankel data for order %f", order));
  endif

  % 1. rescale sinogram data to H-data on [-1,1]
  if (verbose_level >= 1)
    printf("Symmetrization and scaling of data towards H-data...");
  endif
  f_start = tic;
  [h_data, h_grid] = TransformHankelToH(hankel_data, hankel_grid, order, integer_order_flag = false);
  if (verbose_level >= 1)
    printf("Done. \n\tElapsed time in seconds: %.2f\n", toc(f_start));
  endif
  if (do_plotting) % plot transformed data
    figure, ...
    plot(h_grid, h_data), ...
    title(sprintf("H-transformed data for order %f", order));
  endif

  % 2. transform H-data to sinogram
  if (verbose_level >= 1)
    printf("Transforming H-data to a 3D-sinogram...");
  endif
  f_start = tic;
  [sinogram_data, sinogram_theta_grid, sinogram_phi_grid, sinogram_s_grid] = ...
                         reconSinogramFromH3D(h_data, ...
                                              h_grid, ...
                                              theta_grid_size, ...
                                              phi_grid_size, ...
                                              c, ...
                                              npswf, ...
                                              order, ...
                                              rsupp=1.0, ...
                                              verbose_level=verbose_level
                                              );
  if (verbose_level >= 1)
    printf("Done. \n\tElapsed time in seconds: %.2f\n", toc(f_start));
  endif

  % 3. reconstruct 'v' via Radon inversion
  if (verbose_level >= 1)
    printf("Inversion of Radon transform in 3D/reconstruction of 'v'...");
  endif
  f_start = tic;
  s_grid_size = 2*length(hankel_grid)-1; % REFACTOR : dirt when defining this
  [nodes, values, jacobian_weights, Fs] = rtft3d(sinogram_data, ...
                                                 theta_grid_size, ...
                                                 phi_grid_size, ...
                                                 s_grid_size, ...
                                                 1.0, ...
                                                 rpadd, ...
                                                 true, ...
                                                 verbose_level,
                                                 true); ...
  nufft_ngrid = 2*floor(Fs) + 2; % minimal grid size for correct usage of NUFFT
  if (verbose_level >= 1)
    printf("\tNUFFT will reconstruct on image of dimension %dx%d\n", nufft_ngrid, nufft_ngrid);
  endif
  % REFACTOR : normally function must return function and its domain
  v = nfft_reconstruct_3d(nufft_ngrid, nodes, values, jacobian_weights, Fs, verbose_level);
  [v_grid_x, v_grid_y, v_grid_z] = ndgrid(linspace(-1.,1.,nufft_ngrid + 1)(1:nufft_ngrid),
                                          linspace(-1.,1.,nufft_ngrid + 1)(1:nufft_ngrid),
                                          linspace(-1.,1.,nufft_ngrid + 1)(1:nufft_ngrid));
  if (verbose_level >= 1)
    printf("Done. \n\tElapsed time in seconds: %.2f\n", toc(f_start));
  endif
  % plot Radon reconstruction before angle averaging
  if (do_plotting == true)
    mid = floor(nufft_ngrid/2);
    v_slice_Z = squeeze(v(mid, :, :));
    v_slice_Y = squeeze(v(:, mid, :));
    v_slice_X = squeeze(v(:, :, mid));

    figure, imshow(real(v_slice_X), [min(min(real(v_slice_X))), max(max(real(v_slice_X)))]), ...
      colorbar, title(sprintf("middle X-slice, Re(v), c=%f, npswf=%d, ord=%f", c, npswf, order));
    figure, imshow(real(v_slice_Y), [min(min(real(v_slice_Y))), max(max(real(v_slice_Y)))]), ...
      colorbar, title(sprintf("middle Y-slice, Re(v), c=%f, npswf=%d, ord=%f", c, npswf, order));
    figure, imshow(real(v_slice_Z), [min(min(real(v_slice_Z))), max(max(real(v_slice_Z)))]), ...
      colorbar, title(sprintf("middle Z-slice, Re(v), c=%f, npswf=%d, ord=%f", c, npswf, order));

    figure, imshow(imag(v_slice_X), [min(min(imag(v_slice_X))), max(max(imag(v_slice_X)))]), ...
      colorbar, title(sprintf("middle X-slice, Im(v), c=%f, npswf=%d, ord=%f", c, npswf, order));
    figure, imshow(imag(v_slice_Y), [min(min(imag(v_slice_Y))), max(max(imag(v_slice_Y)))]), ...
      colorbar, title(sprintf("middle X-slice, Im(v), c=%f, npswf=%d, ord=%f", c, npswf, order));
    figure, imshow(imag(v_slice_Z), [min(min(imag(v_slice_Z))), max(max(imag(v_slice_Z)))]), ...
      colorbar, title(sprintf("middle X-slice, Im(v), c=%f, npswf=%d, ord=%f", c, npswf, order));
  endif

  % 4. apply radial multiplier 'R' and average over angles on S^2
  if (verbose_level >= 1)
    printf("Angle averaging of reconstructed signal...");
  endif
  if (verbose_level >= 2)
    printf("\n\tnumber of angles for averaging %d\n", phi_grid_size*theta_grid_size);
  endif
  f_start = tic;
  RR = sqrt(v_grid_x.^2 + v_grid_y.^2 + v_grid_z.^2); % normalize to r^1
  [recon, recon_grid] = HarmonicProjection3DToRadial1D( v.*RR, ...
                                                        v_grid_x, ...
                                                        v_grid_y, ...
                                                        v_grid_z,
                                                        recon_grid_size, ...
                                                        n, ...
                                                        phi_grid_size, ...
                                                        theta_grid_size, ...
                                                        interp_method = "linear",
                                                        from_nfft = true
                                                        );
  if (verbose_level >= 1)
    printf("Done. \n\tElapsed time in seconds: %.2f\n", toc(f_start));
  endif

  if (do_plotting == true)
    rrecon = real(recon);
    irecon = imag(recon);

    % plot 1D reconstruction
    figure, plot(recon_grid, rrecon), title(sprintf("Re(f) c=%f, npswf=%d, order=%f", c, npswf, order));
    figure, plot(recon_grid, irecon), title(sprintf("Im(f) c=%f, npswf=%d, order=%f", c, npswf, order));

    % plot 2D embedding
    [XX, YY] = meshgrid(linspace(-1,1,plot2d_size), linspace(-1,1,plot2d_size));
    RR = sqrt(XX.^2 + YY.^2);
    recon2d = zeros(plot2d_size, plot2d_size);
    for ind = 1:plot2d_size % its ugly but fast
      recon2d(ind, :) = interp1(recon_grid, recon, RR(ind, :), "linear", 0.0);
    endfor
    rrecon2d = real(recon2d);
    irecon2d = imag(recon2d);
    figure, imshow(rrecon2d, [min(min(rrecon2d)), max(max(rrecon2d))]), title(sprintf("2D Re(f) c=%f, npswf=%d, order=%f", c, npswf, order)), colorbar;
    figure, imshow(irecon2d, [min(min(irecon2d)), max(max(irecon2d))]), title(sprintf("2D Im(f) c=%f, npswf=%d, order=%f", c, npswf, order)), colorbar;
  endif
  if (verbose_level >= 1)
    printf("Reconstruction finished. Total elapsed time in seconds: %.2f\n", toc(start));
  endif

endfunction

