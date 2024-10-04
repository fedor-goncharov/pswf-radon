% reconstruction from Hankel transform of integer order limited on [0, c]
% via PSWFs SVD decomposition and Radon inversion in 2D

function [recon, recon_grid] = reconFromHankelPSWF2D(hankel_data, ...
                                                 hankel_grid, ...
                                                 c, ...
                                                 order,
                                                 npswf, ...
                                                 phi_grid_size,
                                                 recon_grid_size,
                                                 interp_method = "linear", ...
                                                 rpadd = 1.3,
                                                 verbose_level = 0, ...
                                                 do_plotting = false, ...
                                                 plot2d_size = 64
                                                 )
  % depends NFFT library (see Chemnitz-TU for NFFT)
  % Script performs reconstruction of a function from its Fourier transforms at 'nodes'
  % with values 'values'. Mathematically it works as a discretized version of (inverse) Fourier integral
  % in 3D. Result is given as a 3D matrix of size : ngrid x ngrid x ngrid, which, in turn, is a
  % grid on [-1,1)x[-1,1)x[-1,1).

  % Usage of the script
  % ngrid : number of points (as intervals) in [-1.0, 1.0)
  % nodes : matrix of size (number_of_nodes x 3), where '3' stands for [x, y, z], these are the
  %         points in space where Fourier transform is known
  % values: vector of size (number_of_nodes) values of Fourier transform of a function in 'nodes'
  % jacobian_weights : volumes at 'nodes' in the Riemann summ of discretized Fourier integral
  % Fs    : double of bandwidth of data (for NFFT usage it is necessary 'ngrid > 2*Fs')

  % Description:
  %   ...
  % params:
  %   ...
  % assumed:
  %   reconstructed signal is supported on [0, 1]

  % input checks
  assert(hankel_grid(1) == 0.0 && hankel_grid(end) == c, "Problem with grid - check doc for correct usage. Abort. \n");
  assert(isinteger(npswf) && npswf >= uint8(0), "Provided 'npswf' is not of integer type or negative (e.g. npswf=uint8(X)). Abort.\n");
  assert(isinteger(order), "Provided order is not of integer type. Did you cast it to uintX(Y)?\n");
  assert(order >= uint8(0), "Provided order is negative. Abort.\n")

  order = double(order);

  % global timer
  start = tic;

  if (verbose_level >= 1)
    printf("Reconstruction from Hankel transform data of signal supported in [0,1]\n");
    printf("Main parameters\n");
    printf("\tData length: %d, grid start: %f, grid end: %f\n", length(hankel_data), hankel_grid(1), hankel_grid(end));
    printf("\tBandwidth 'c': %f\n", c);
    printf("\tHankel transform order: %d\n", order);
    printf("\tNumber of PSWFs used: %d\n", npswf);
    printf("\tNumber of projections for Radon embedding: %d\n", phi_grid_size);
    printf("\tPad Radon data with zeros up to: %f\n", rpadd);
    printf("Auxiliary parameters\n");
    printf("\tInterpolation method for angle averaging: %s\n", interp_method);
    printf("\tDo plots: %d\n", do_plotting);
  endif

  if (do_plotting) % plot data
    figure, ...
    plot(hankel_grid, hankel_data), ...
    title(sprintf("Hankel data for order %d", order));
  endif

  % 1. rescale data to [-1,1] to H-data
  if (verbose_level >= 1)
    printf("Symmetrization and scaling of data...");
  endif
  f_start = tic;
  [h_data, h_grid] = TransformHankelToH(hankel_data, hankel_grid, order, integer_order_flag = true);
  if (verbose_level >= 1)
    printf("Done. \n\tElapsed time in seconds: %.2f\n", toc(f_start));
  endif
  if (do_plotting) % plot transformed data
    figure, ...
    plot(h_grid, h_data), ...
    title(sprintf("H-transformed data for order %d", order));
  endif

  % 2. reconstruct sinogram from H-data
  if (verbose_level >= 1)
    printf("Reconstruction of sinogram data from H-data...");
  endif
  f_start = tic;
  [sinogram_data, ...
   sinogram_phi_grid, ...
   sinogram_s_grid] = reconSinogramFromH2D( h_data, ...
                                            h_grid, ...
                                            phi_grid_size, ...
                                            c, ...
                                            npswf, ...
                                            order, ...
                                            1.0, ...
                                            verbose_level=verbose_level
                                            );
  if (verbose_level >= 1)
    printf("Done. \n\tElapsed time in seconds: %.2f\n", toc(f_start));
  endif

  if (do_plotting == true)
    rsino = real(sinogram_data);
    isino = imag(sinogram_data);
    figure, imshow(rsino, [min(min(rsino)), max(max(rsino))]), ...
      title(sprintf("Re(sinogram) order=%d, c=%.2f, nphi=%d", order, c, phi_grid_size));
    figure, imshow(isino, [min(min(isino)), max(max(isino))]), ...
      title(sprintf("Im(sinogram) order=%d, c=%.2f, nphi=%d", order, c, phi_grid_size));
  endif

  % 4. reconstruct 'v' via Radon inversion
  if (verbose_level >= 1)
    printf("Inversion of Radon transform in 2D/reconstruction of 'v'...");
  endif
  f_start = tic;
  [v_nodes, v_values, v_jacobian_weights, Fs] = rtft2d(sinogram_data, ...
                                                 size(sinogram_data, 1), ...
                                                 size(sinogram_data, 2), ...
                                                 1.0, ...
                                                 rpadd-1.0, ...
                                                 fourier_padding_flag = true, ...
                                                 verbose_level = verbose_level
                                                );
  % NUFFT
  nufft_ngrid = 2*floor(Fs) + 2; % minimal grid size for correct usage of NUFFT
  v_data = nfft_reconstruct_2d(nufft_ngrid, v_nodes, v_values, v_jacobian_weights, Fs, verbose_level=verbose_level);
  [v_grid_y, v_grid_x] = ndgrid(linspace(-1.,1.,nufft_ngrid + 1)(1:nufft_ngrid),
                                linspace(-1.,1.,nufft_ngrid + 1)(1:nufft_ngrid));

  % IRADON (slower) if to activae - change: 'to2pi=false' in 'TransformHToFcSinogram2D', 'grid_pos_o = false' in S^1 averaging
  %padd_to_one_size = floor((sqrt(2)-1.0)*size(sinogram_data)(2)/2) + 1;
  %padd_to_extra_size = (rpadd-1.0 > 0.0)*floor((rpadd - 1.0)*(size(sinogram_data)(2)/2 + padd_to_one_size)) + 1;
  %if (verbose_level >= 2)
  %  printf("\n\tzero padding of size %d on each size.\n", padd_to_one_size + padd_to_extra_size);
  %endif
  %[v_data, v_grid_x, v_grid_y] = TransformSinogramToV2D(sinogram_data, ...
  %                                                      sinogram_s_grid, ...
  %                                                      sinogram_phi_grid, ...
  %                                                      padd_to_one_size + padd_to_extra_size);
  if (verbose_level >= 1)
    printf("Done. \n\tElapsed time in seconds: %.2f\n", toc(f_start));
  endif

  % plot Radon reconstruction before angle averaging
  if (do_plotting == true)
    figure, imshow(real(v_data), [min(min(real(v_data))), max(max(real(v_data)))]), colorbar, title(sprintf("Re(v) c=%f, npswf=%d, order=%f", c, npswf, order));
    figure, imshow(imag(v_data), [min(min(imag(v_data))), max(max(imag(v_data)))]), colorbar, title(sprintf("Im(v) c=%f, npswf=%d, order=%f", c, npswf, order));
    v_data_dimsize = size(v_data)(1);
    extended_limit = v_grid_x(end, end);
    figure, plot(linspace(-extended_limit, extended_limit, v_data_dimsize), real(v_data(ceil(v_data_dimsize/2), :))), ...
      title(sprintf("central section Re(v) c=%f, npswf=%d, order=%f", c, npswf, order));
    figure, plot(linspace(-extended_limit, extended_limit, v_data_dimsize), imag(v_data(ceil(v_data_dimsize/2), :))), ...
      title(sprintf("central section Im(v) c=%f, npswf=%d, order=%f", c, npswf, order));
  endif

  % 4. apply radial multiplier sqrt(R) and average over angles
  if (verbose_level >= 1)
    printf("Angle averaging of reconstructed signal...");
  endif
  if (verbose_level >= 2)
    printf("\n\tnumber of angles for averaging %d\n", phi_grid_size);
  endif
  f_start = tic;
  sqRR = sqrt(sqrt(v_grid_x.^2 + v_grid_y.^2)); % normalize to the square root % DEBUG
  [recon, recon_grid] = HarmonicProjection2DToRadial1D(v_data.*sqRR, ...
                                                       v_grid_x, ...
                                                       v_grid_y, ...
                                                       recon_grid_size, ...
                                                       order, ...
                                                       phi_grid_size, ...
                                                       interp_method = interp_method, ...
                                                       divide_2pi = true, ...
                                                       grid_pos_o = true);
  if (verbose_level >= 1)
    printf("Done. \n\tElapsed time in seconds: %.2f\n", toc(f_start));
  endif

  if (do_plotting == true)
    rrecon = real(recon);
    irecon = imag(recon);

    % plot 1D
    figure, plot(recon_grid, rrecon), title(sprintf("Re(f) c=%f, npswf=%d, order=%f", c, npswf, order));
    figure, plot(recon_grid, irecon), title(sprintf("Im(f) c=%f, npswf=%d, order=%f", c, npswf, order));

    % plot 2D
    [XX, YY] = meshgrid(linspace(-1,1,plot2d_size), linspace(-1,1,plot2d_size));
    RR = sqrt(XX.^2 + YY.^2);
    recon2d = zeros(plot2d_size, plot2d_size);
    for ind = 1:plot2d_size % ugly but fast
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

