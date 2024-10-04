% rtft2d.m
% function evaluates numerically Fourier integral of the sinogram
% in order to obtain Fourier data in 2D; mapping between Fourier transform
% of the sinogram and original signal is described by 'Projection Theorem'

function [nodes, values, jacobian_weights, Fs] = rtft2d(sinogram, ...
                                                    nphi, ...
                                                    nshift, ...
                                                    rsupp, ...
                                                    padding = 1.0, ...
                                                    fourier_padding_flag = true, ...
                                                    verbose_level = 1, ...
                                                    apply_cosine_filter = true
                                                   )
% depends on cartprod.m

% This script reads data given by ray transforms in 2D from file and performs
% 1D Fourier transform along shift argument. Returns values (values) of the
% 1D Fourier integral and respective nodes (nodes) with volumes corresponding to nodes (jacobian_weights).

% Usage of the script:

%   nphi              : number of projections in azimuth angle
%   nshift            : number of rays per one direction, shift [-1, 1]
%   rsupp             : radius of the support of the test-function
%   padding           : multiplying factor for appending data with zeros

% It is assumed that angles 'phi' and 'shifts' are uniformly spaced
% in their intervals (0,2pi], [-1,1], respectively.

% OUTPUT : nodes - nodes in 3D spaces, where Fourier Transform of the function is known
%          values - values of FT of the function
%          jacobian_weights - volume of cells for each node in fourier space
%          Fs - twice maximal frequency of the grid to which nodes belong (useful for further NFFT reconstructions)



  % set grid in Radon space
  phi_vec = vec(linspace(0, 2*pi, nphi + 1)(1:end-1));     % directions 'phi' on circle [0, 2*pi)

  % set spatial and frequency grid sizes
  dphi   = 2.*pi/nphi;
  dshift = 2.*rsupp / (nshift-1);             % time step of the signal
  padding_size = floor(padding * nshift);     % zero padding of projections ('padding_size' added on each side of projection, i.e., 2 times)
  ntotal = nshift + 2*padding_size;           % total points for given projection;
  Fs = 1./dshift;                             % by DFT representation of Fourier integral full bandwidth of projection is [0, Fs)
  dfreq  = Fs/ntotal;                         % discretization step in frequency domain

  % set frequencies of 1d projections
  frequencies_centered = vec(linspace(-Fs/2, Fs/2, ntotal + 1)(1:end-1)); % endpoint = false

  % set output containers
  per_phi_vsize = prod(size(frequencies_centered));  % size of a Fourier slice for given 'phi'

  nodes = zeros(per_phi_vsize*nphi, 2);                                         % 3d frequencies
  values = zeros(per_phi_vsize*nphi, 1);                                        % values of Fourier integrand in nodes
  jacobian_weights = zeros(per_phi_vsize*nphi, 1);                              % jacobian multipliers

  if verbose_level >= 1
    printf("rtft2d :: Generating nodes in Fourier space, computing jacobians and integrand values...\n");
  endif

  for i_phi = 1 : nphi
      if verbose_level >= 2
        printf("\titeraion %d of %d", i_phi, nphi);
      endif
      % append nodes
      phi = phi_vec(i_phi);
      direction = [cos(phi) sin(phi)];
      nodes((i_phi-1)*per_phi_vsize + 1:i_phi*per_phi_vsize, :) = frequencies_centered*direction;

      % append jacobian
      jacobian_weight = 0.5*dfreq*dphi*abs(frequencies_centered);

      if (apply_cosine_filter == true)
      %raised_cosine_filter = 0.5*(1.0 + cos(pi*frequencies_centered/Fs));
        raised_cosine_filter = sinc(2*frequencies_centered/Fs);
        jacobian_weight = jacobian_weight .* raised_cosine_filter;
      endif

      % append jacobian weights of cells
      jacobian_weights((i_phi-1)*per_phi_vsize + 1:i_phi*per_phi_vsize) = jacobian_weight;

      % 1D Fourier integral via 1D FFT
      projection = vec(sinogram(i_phi, :));
      fft_vec = fftshift(fft(ifftshift(padarray(projection, padding_size))));
      fft_vec = dshift*fft_vec;
      values((i_phi-1)*per_phi_vsize + 1:i_phi*per_phi_vsize) = fft_vec;
  end

  if verbose_level >= 1
    printf("rtft2d :: Done. %d primary nodes created.\n", size(values, 1));
  endif

  % stabilization - append nodes between ball of radius Fs/2 and [-Fs/2,Fs/2]^2
  size_append = 0;
  if fourier_padding_flag == true
    if verbose_level >= 1
      printf("rtft2d :: Appending exterior nodes...");
    endif
    % compute nodes
    append_nodes = cartprod(frequencies_centered, frequencies_centered);
    append_nodes = [append_nodes, sqrt(sum(append_nodes.^2, 2))];
    append_nodes = append_nodes(append_nodes(:, 3) > Fs/2, :);
    append_nodes = append_nodes(:, [1 2]);
    % append
    nodes = [nodes; append_nodes];
    size_append = size(append_nodes, 1);
    if verbose_level >= 1
      printf("rtft2d :: Done. %d nodes were appended.\n", size_append);
    endif
  endif

  % append extra zeros to values and 'uniform grid' jacobians if necessary
  values = [values; zeros(size_append, 1)];
  jacobian_weights = [jacobian_weights; ones(size_append, 1)*(dfreq^2)];

  if verbose_level >= 1
    printf("rtft2d :: All done.\n");
  endif
end

