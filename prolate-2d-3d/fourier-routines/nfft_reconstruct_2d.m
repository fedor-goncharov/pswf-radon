% Reconstruction of a test function in 2D from its Fouirer transforms

function rec = nfft_reconstruct_2d(ngrid, nodes, values, jacobian_weights, Fs, verbose_level=1, num_threads=-1)
% depends NFFT library from TU Chemnitz

% Script performs reconstruction of a function from its Fourier transforms at 'nodes'
% with values 'values'. Reconstruction is performed via discretized inverse Fourier transform
% in 2D. Result is given as a 2D matrix of size : ngrid x ngrid, which, in turn, is a
% grid on [-1,1)x[-1,1).

% Usage
% ngrid : number of points (as intervals) in [-1.0, 1.0)
% nodes : matrix of size (number_of_nodes x 2), where '2' stands for [x, y], these are the
%          points in space where Fourier transform is known
% values : vector of size (number_of_nodes) values of Fourier transform of a function in 'nodes'
% jacobian_weights : volumes at 'nodes' in the Riemann summ of discretized Fourier integral
% Fs    : double of bandwidth of data (for NFFT usage it is necessary that 'ngrid > 2*Fs')

  if (Fs/ngrid >= 0.5)
    error(sprintf("must have Fs/ngrid < 0.5 but recieved %f => nodes domain cannot be larger than NFFT domain [-1/2, 1/2)^3", Fs/ngrid));
  endif

  % jacobian weightening
  summands = jacobian_weights .* values;


  %normalization of 'nodes' so that they belong to [-1/2, 1/2)^3 (see theory)
  if (verbose_level >= 1)
    printf("nfft_reconstruct_2d :: scaling of nodes to comply with NUFFT routines...");
  endif
  normalizer = 2.0/ngrid;
  nodes_normalized = nodes*normalizer;
  nodes_normalized = nodes_normalized';
  if (verbose_level >= 1)
    printf("Done.\n");
  endif

  if (verbose_level >= 1)
    printf("nfft_reconstruct_2d :: NUFFT initializations (plan, data, nodes, etc.)...");
  endif
  if (num_threads > 0) % if num_threads < 0 use maximum number of threads
    nfft_set_num_threads(num_threads);
  endif
  n=2^(ceil(log(ngrid)/log(2))+1);
  plan = nfft_init_guru(2, ngrid, ngrid, size(nodes_normalized,2), n, n, 8,
  NFFT_OMP_BLOCKWISE_ADJOINT, FFTW_ESTIMATE);

  % set nodes in plan
  nfft_set_x(plan, nodes_normalized);

  % precomputations
  nfft_precompute_psi(plan);

  % set Fourier coefficients
  nfft_set_f(plan, summands);

  if (verbose_level >= 1)
    printf("nfft_reconstruct_2d :: Done.\n");
  endif

  % inversion
  if (verbose_level >= 1)
    printf("nfft_reconstruct_2d :: Inversion via adjoint transform...");
  endif
  nfft_adjoint(plan);
  if (verbose_level >= 1)
    printf("Done.\n");
  endif

  % return
  rec = reshape(nfft_get_f_hat(plan), ngrid, ngrid);

  % clean memory
  if (verbose_level >= 1)
    printf("nfft_reconstruct_2d :: NUFFT finalization ...");
  endif
  nfft_finalize(plan);
  if (verbose_level >= 1)
    printf("Done.\n");
  endif
end
