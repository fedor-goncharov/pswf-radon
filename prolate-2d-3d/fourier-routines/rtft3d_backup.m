% rtft3d.m

function [nodes, values, jacobian_weights, Fs] = rtft3d(rt_matrix, ...
                                                    ntheta, ...
                                                    nphi, ...
                                                    nshift, ...
                                                    rsupp, ...
                                                    padding, ...
                                                    fourier_padding_flag = true)
% depends on lgwt.m, cartprod.m. centered_frequencies.m

% This script reads 3D-Radon transforms (over 2D planes) from file and computes
% 1D Fourier transforms along all projections. The returned result consists of nodes of the
% aforementioned 1D Fourier integral, its values and
% volume size (jacobian weights) of each node.

% Usage of the script

% filename : file where the data is stored in the form [sigma, phi, theta]
% ntheta   : number of projections in polar angle (0, pi)
% nphi     : number of projections in azimuth angle [0, 2*pi)
% nshift   : number of hyperplanes per one direction
% rsupp    : radius of the support, where the test function is lying
% padding_coeff : multiplcative factor for zero padding for each projection for FFT's
%                 (by default=4, i.e., 4 times the length of shift array)

% Geometry

% uniform distribution of 'phi', Gaussian polar angles 'theta'
% (i.e. theta_j = arccos(t_j), (t_j, j = 1, ntheta) - Gaussian points on [-1, 1])
% and uniform 'shifts' on [-1,1]

% OUTPUT : nodes - nodes in 3D spaces, where FT of the function is known
%          values - values of FT of the function
%          jacobian_weights - volume of cells for each node in fourier space


  % set grid in Radon space
  phi_vec = vec(linspace(0, 2*pi, nphi + 1)(1:end-1));     % directions 'phi' on circle [0, 2*pi)
  u_vec = vec(linspace(-1,1,ntheta + 2)(2:end-1));
  theta_vec = acos(u_vec);                                 % directions 'theta' on (0, pi)

  % set spatial and frequency grid sizes
  du = u_vec(2) - u_vec(1);
  dphi   = 2.*pi / nphi;
  dshift = 2.*rsupp / (nshift-1);             % time step of the signal
  padding_size = floor(padding * nshift);     % zero padding of projections ('padding_size' added on each side of projection, i.e., 2 times)
  ntotal = nshift + 2*padding_size;           % total points for given projection;

  Fs = 1./dshift                              % by DFT representation of Fourier integral full bandwidth of projection is [0, Fs)
  dfreq  = Fs/ntotal;                         % discretization step in frequency domain

  % set frequencies of 1d projections
  frequencies_centered = vec(linspace(-Fs/2, Fs/2, ntotal + 1)(1:end-1)); % endpoint = false

  % DEBUG
  %frequencies = vec(0:ntotal-1)*dfreq;
  %frequencies_centered = vec(centered_frequencies(ntotal, dfreq)); % idk why a special function is needed in this case
  %length(frequencies_centered)/2
  %max(frequencies_centered)

  % generate Nyquist sinc filter
  %half_nyq_frequency = ntotal*dfreq / 2;
  %sinc_filter = (sinc(frequencies_centered / half_nyq_frequency)).^2;

  % set containers
  per_theta_vsize = prod(size(frequencies_centered))*nphi;                      % size of a Fourier slice for given 'theta'
  nodes = zeros(ntheta*per_theta_vsize, 3);                                     % array of nodes in frequency domain
  values = zeros(ntheta*per_theta_vsize, 1);                                    % array of values of Fourier integral in nodes
  jacobian_weights = zeros(ntheta*per_theta_vsize, 1);                          % array of jacobian multipliers for each node

  % -------------------- Primary nodes in Fourier grid -------------------------
  printf("Generating nodes in Fourier space...\n");

  for i_theta = 1 : ntheta
    printf("\tgeneration of Fourier nodes : iteration (of %d) : %d\n", ntheta, i_theta);

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
  printf("Done.\n");
  % ----------------------------------------------------------------------------

  % ------------------Appending secondary nodes in Fourier space ---------------
  % ----------------- needed to stabilize gridding of NFFT in future -----------
  size_append = 0;
  if (fourier_padding_flag == true)
    printf("Primary nodes created. Appending secondary nodes for NFFT stabilization...\n");
    % append zero nodes outside the ball of radius Fs/2 up to cube [-Fs/2, Fs/2)^3
    append_nodes = cartprod(frequencies_centered, ...
                            frequencies_centered, ...
                            frequencies_centered);                 % (length(frequencies_centered)^3, 3)
    append_nodes = [append_nodes, sqrt(sum(append_nodes.^2, 2))];  % 3d-vector norm in last column
    append_nodes = append_nodes(append_nodes(:, 4) > Fs/2, :);     % choose nodes outside the ball
    append_nodes = append_nodes(:, [1 2 3]);                       % forget norm
    nodes = [nodes; append_nodes];

    size_append = size(append_nodes, 1);  % track to add extra jacobian, value rows
    printf("Done. %d secondary nodes added\n", size_append);
  endif
  % ----------------------------------------------------------------------------

  % ---------------------- Computation of jacobian weights ---------------------
  printf("Computing jacobian weights for all nodes...");
  for i_theta = 1 : ntheta
    printf("\tjacobian weights : iteration (of %d) : %d\n", ntheta, i_theta);
    jacobian_weight = 0.5 * dfreq * dphi * (frequencies_centered.^2) * du;

    % add smoothing filter
    % TODO
    % legacy jacobian_weight = jacobian_weight .* sinc_filter;

    jacobian_weight = repmat(jacobian_weight, nphi, 1);
    jacobian_weights(((i_theta-1)*per_theta_vsize + 1):i_theta*per_theta_vsize) = jacobian_weight;
  endfor

  % add 'ones' for jacobian weights at nodes out of the ball
  jacobian_weights = [jacobian_weights; ones(size_append, 1)*(dfreq^3)];
  printf("Done.\n");
  % ----------------------------------------------------------------------------


  % ---------Evaluation of 1d-Fourier transform for each projection ------------
  printf("Computing values of for the Fourier transform at nodes...");

  for i_theta = 1 : ntheta
    printf("\tFourier transforms at nodes : iteration (of %d) : %d\n", ntheta, i_theta);
    % fft along dimension 'nshift'
    rt_theta_slice =  permute(rt_matrix(i_theta, :, :), [3 2 1]);                        % ntotal x nphi
    fft_theta_slice = fft(ifftshift(padarray(rt_theta_slice, [padding_size 0]), 1));     % fft along vertical dimension
    fft_theta_slice = dshift * fftshift(fft_theta_slice , 1);                            % fftshift along vertical dimension

    % reshaping slice to vector
    values_add = reshape(fft_theta_slice, nphi*ntotal, 1);
    values(((i_theta-1)*per_theta_vsize + 1):i_theta*per_theta_vsize) = values_add;
  endfor
  % add zeros at nodes out of the ball
  values = [values; zeros(size_append, 1)];
  printf("Done.\n");
  %-----------------------------------------------------------------------------

endfunction
