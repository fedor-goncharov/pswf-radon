% rtft3d.m
% function evaluates numerically Fourier integral of the sinogram
% in order to obtain Fourier data in 3D; mapping between Fourier transform
% of the sinogram and original signal is described by 'Projection Theorem'

function [nodes, values, jacobian_weights, Fs] = rtft3d(sinogram, ...
                                                    ntheta, ...
                                                    nphi, ...
                                                    nshift, ...
                                                    rsupp, ...
                                                    padding, ...
                                                    fourier_padding_flag = true, ...
                                                    verbose_level = 1, ...
                                                    apply_cosine_filter=true
                                                    )
% depends on cartprod.m (Matlab File Exchange)

% This script reads 3D-Radon transforms (over 2D planes) from file and computes
% 1D Fourier transforms along all projections. The returned result consists of nodes of the
% aforementioned 1D Fourier integral, its values and volume sizes (jacobian weights).

% Usage of the script

% sinogram : matrix of values of shape (theta, phi, shift)
% ntheta   : number of projections in polar angle (0, pi), arccos-uniform, no endpoints
% nphi     : number of projections in azimuth angle [0, 2*pi), uniform, no endpoint
% nshift   : number of hyperplanes per projection
% rsupp    : radius of the support, where the test function is lying
% padding  : multiplcative factor for zero padding for each projection for FFT's
% fourier_padding_flag : TODO ...
% verbose_level : TODO ...

% OUTPUT : nodes - nodes in 3D spaces, where FT of the function is known
%          values - values of FT of the function
%          jacobian_weights - volume of cells for each node in fourier space
%          Fs - twice maximal frequency of the grid to which nodes belong (useful for further NFFT reconstructions)


  % set grid in Radon space
  phi_vec = vec(linspace(0, 2*pi, nphi + 1)(1:end-1));     % directions 'phi' on circle [0, 2*pi)
  u_vec = vec(linspace(-1,1,ntheta + 2)(2:end-1));
  theta_vec = acos(u_vec);                                 % directions 'theta' on (0, pi)
  du = u_vec(2) - u_vec(1);
  dphi   = 2.*pi / nphi;
  dshift = 2.*rsupp / (nshift-1);             % time step of the signal
  padding_size = floor(padding * nshift);     % zero padding of projections ('padding_size' added on each side of projection, i.e., 2 times)
  ntotal = nshift + 2*padding_size;           % total points for given projection;

  Fs = 1./dshift;                             % by DFT representation of Fourier integral full bandwidth of projection is [0, Fs)
  dfreq  = Fs/ntotal;                         % discretization step in frequency domain

  % set frequencies of 1d projections
  frequencies_centered = vec(linspace(-Fs/2, Fs/2, ntotal + 1)(1:end-1)); % endpoint = false

  % set output containers
  per_theta_vsize = prod(size(frequencies_centered))*nphi;                      % size of a Fourier slice for given 'theta'
  nodes = zeros(ntheta*per_theta_vsize, 3);                                     % array of nodes in frequency domain
  values = zeros(ntheta*per_theta_vsize, 1);                                    % array of values of Fourier integral in nodes
  jacobian_weights = zeros(ntheta*per_theta_vsize, 1);                          % array of jacobian multipliers for each node

  % -------------------- Primary nodes in Fourier grid -------------------------
  if verbose_level >= 1
    printf("rtft3d :: Generating nodes in Fourier space...\n");
  endif

  for i_theta = 1 : ntheta
    if verbose_level >= 2
      printf("\tgeneration of Fourier nodes : iteration (of %d) : %d\n", ntheta, i_theta);
    endif

    % make of directions for given theta
    theta = theta_vec(i_theta);
    direction_theta = reshape([sin(theta)*cos(phi_vec)'; ...
                               sin(theta)*sin(phi_vec)'; ...
                               cos(theta)*ones(1,nphi)], 1, []);

    % 'direction' * frequency makes a node
    % 1. set nodes in matrix (frequencies, (dir1.x, dir1.y, dir1.z, dir2.x, ...))
    % 2. reshape nodes into matrix [freq1.x, freq1.y, freq1.z;
    %                               freq2.x, freq2.y, freq2.z]
    %                              of size (frequencies*nphi, 3)
    nodes_theta = frequencies_centered * direction_theta;
    nodes_to_add = permute(reshape(nodes_theta, length(frequencies_centered), 3, nphi), [1 3 2]);

    % save to container
    nodes_to_add = reshape(nodes_to_add, [], 3);
    nodes(((i_theta-1)*per_theta_vsize + 1):i_theta*per_theta_vsize, :) = nodes_to_add;
  endfor
  if verbose_level >= 1
    printf("rtft3d :: Done.\n");
  endif
  % ----------------------------------------------------------------------------

  % ------------------Appending secondary nodes in Fourier space ---------------
  % ----------------- needed to stabilize gridding of NFFT in future -----------
  size_append = 0;
  if (fourier_padding_flag == true)
    if verbose_level >= 1
      printf("rtft3d :: Primary nodes created. Appending secondary nodes for NFFT stabilization...\n");
    endif
    % append zero nodes outside the ball of radius Fs/2 up to cube [-Fs/2, Fs/2)^3
    append_nodes = cartprod(frequencies_centered, ...
                            frequencies_centered, ...
                            frequencies_centered);                 % (length(frequencies_centered)^3, 3)
    append_nodes = [append_nodes, sqrt(sum(append_nodes.^2, 2))];  % 3d-vector norm in last column
    append_nodes = append_nodes(append_nodes(:, 4) > Fs/2, :);     % choose nodes outside the ball
    append_nodes = append_nodes(:, [1 2 3]);                       % forget norm
    nodes = [nodes; append_nodes];

    size_append = size(append_nodes, 1);  % track to add extra jacobian, value rows
    if verbose_level >= 1
      printf("rtft3d :: Done: %d secondary nodes added.\n", size_append);
    endif
  endif
  % ----------------------------------------------------------------------------

  % ---------------------- Computation of jacobian weights ---------------------
  if verbose_level >= 1
    printf("rtft3d :: Computing jacobian weights for all nodes...");
  endif
  for i_theta = 1 : ntheta
    if verbose_level >= 2
      printf("\tjacobian weights : iteration (of %d) : %d\n", ntheta, i_theta);
    endif
    jacobian_weight = 0.5 * dfreq * dphi * (frequencies_centered.^2) * du;

    % raised cosine smoothing filter in Fourier domain
    if (apply_cosine_filter == true)
      %raised_cosine_filter = 0.5*(1.0 + cos(pi*frequencies_centered/Fs));
      raised_cosine_filter = sinc(2*frequencies_centered/Fs);
      jacobian_weight = jacobian_weight .* raised_cosine_filter;
    endif

    jacobian_weight = repmat(jacobian_weight, nphi, 1);
    jacobian_weights(((i_theta-1)*per_theta_vsize + 1):i_theta*per_theta_vsize) = jacobian_weight;
  endfor

  % add 'ones' for jacobian weights at nodes out of the ball
  jacobian_weights = [jacobian_weights; ones(size_append, 1)*(dfreq^3)]; % TODO : optimize
  if verbose_level >= 1
    printf("rtft3d :: Done.\n");
  endif
  % ----------------------------------------------------------------------------

  % ---------Evaluation of 1d-Fourier transform for each projection ------------
  if verbose_level >= 1
    printf("rtft3d :: Computing values of for the Fourier transform at nodes...");
  endif

## ------------------------------ DEBUG ----------------------------------------
##  for i_theta = 1 : ntheta
##    if verbose_level >= 2
##      printf("\tFourier transforms at nodes : iteration (of %d) : %d\n", ntheta, i_theta);
##    endif
##    % fft along dimension 'nshift'
##    rt_theta_slice =  permute(sinogram(i_theta, :, :), [3 2 1]);                         % result has size (ntotal, nphi)
##    fft_theta_slice = fft(ifftshift(padarray(rt_theta_slice, [padding_size 0]), 1));     % fft along vertical dimension (see doc)
##    fft_theta_slice = dshift * fftshift(fft_theta_slice , 1);                            % fftshift along vertical dimension
##
##    % reshaping slice to vector
##    values_add = reshape(fft_theta_slice, nphi*ntotal, 1);
##    values(((i_theta-1)*per_theta_vsize + 1):i_theta*per_theta_vsize) = values_add;
##  endfor
## ------------------------------- DEBUG ---------------------------------------

  ft_sinogram = fftshift(fft(ifftshift(padarray(sinogram, [0 0 padding_size]), ...
                             3), ...
                         [], 3), ...
                3)*dshift;
  values = reshape(permute(ft_sinogram, [3 2 1]), [], 1); % think of it ;)

  % add zeros at nodes out of the ball
  values = [values; zeros(size_append, 1)];
  if verbose_level >= 1
    printf("rtft3d :: Done.\n");
  endif
  %-----------------------------------------------------------------------------

endfunction
