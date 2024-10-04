% 1d inverse of Fc via SVD

function g = invFc(f, f_grid, c, npswf, verbose=false, pswf_data_path='saved')
   % Description:
   %   inversion of Fc in dimension 1 via SVD series
   % params:
   %   npswf : number of PSWFs used
   %   c     : bandwidth parameter for PSWFs
   % output:
   %   approximate inverse transform of 'f' on grid 'f_grid'
   % assumed:
   %   grid of pswfs is on [-1,1] with endpoints
   %   'f_grid' is on [-1,1] with endpoints and is coarser of pswfs

   g = zeros(1, length(f_grid));
   if verbose
     printf("invFc -> Number of pswfs : %d\n", npswf);
     printf("invFc -> Paramter c : %2.2f\n", c);
   endif
   for n=0:npswf
     % get pswf function
     pswf_struct = pro_open_everything(pswf_data_path, c, 0, n);
     pswf_func = double(pswf_struct.S1.S1); % DOUBLE CAST
     pswf_grid = linspace(-1,1, length(pswf_func)); % uniform grid
     pswf_fc_eigenval = pswf_struct.coefficients.mu;

     % compute dot product between f and pswf
     if (length(f_grid) > length(pswf_grid))
       error(sprintf("invFc -> grid of integrand is finer than of pswf function of order %d\n", n));
     endif
     f_interp = interp1(f_grid, f, pswf_grid, "linear"); % cast function to the grid of pswfs
     dx = pswf_grid(2) - pswf_grid(1);
     l2projection = IntegrateSimpson1D(f_interp.*pswf_func, dx);

     % add new reconstruction term
     g += l2projection/pswf_fc_eigenval*interp1(pswf_grid, pswf_func, f_grid, "linear");
     if verbose
       printf("\tinvFc -> pswf %d -> eigenvalue : %s\n", n, num2str(pswf_fc_eigenval));
       printf("\tinvFc -> pswf %d -> l2-norm,  max, min, std : %f, %f, %f, %f\n", n, dx*norm(g), max(g), min(g), std(g));
     endif
   endfor
endfunction
