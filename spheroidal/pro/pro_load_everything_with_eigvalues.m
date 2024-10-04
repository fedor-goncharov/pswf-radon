%
% Copyright (c) 2014, Ross Adelman, Nail A. Gumerov, and Ramani Duraiswami
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%

%
% Unzip all of the ZIP'd-up prolate spheroidal coefficients and wave functions,
% and organize them into a nice, easily-accessible MATLAB struct.  This is
% where the best combination of prolate spheroidal radial functions is selected
% to minimize the error in the Wronskian.
%
% Arguments:
%     path - the directory in which the prolate spheroidal coefficients and
%            wave functions have been saved and ZIP'd up
%     c - k * a, where k is the wavenumber and 2a is the interfocal distance
%     m - one modal value
%     n - the other modal value, which is equal to m, m + 1, ...
%     name - the same name from pro_calculate_functions
% Return Values:
%     everything - a MATLAB struct with all of the prolate spheroidal
%                  coefficients and wave functions
%
function everything = pro_load_everything_with_eigvalues(path, saved, c, m, n, name, use_vpa = false, vpa_digits = 200)

	everything = struct();
	everything.c = c;
	everything.m = m;
	everything.n = n;
	try
    everything.coefficients = struct();
    unzip(sprintf('%s/pro_%s.zip', path, generate_name(c, m, n)), path);
    if (use_vpa == true)
      everything.coefficients.lambda = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_lambda.txt', path, generate_name(c, m, n)), vpa_digits);
      everything.coefficients.log_abs_lambda = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_log_abs_lambda.txt', path, generate_name(c, m, n)), vpa_digits);
      everything.coefficients.dr = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_dr.txt', path, generate_name(c, m, n)), vpa_digits).';
      everything.coefficients.log_abs_dr = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_log_abs_dr.txt', path, generate_name(c, m, n)), vpa_digits).';
      everything.coefficients.dr_neg = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_dr_neg.txt', path, generate_name(c, m, n)), vpa_digits).';
      everything.coefficients.log_abs_dr_neg = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_log_abs_dr_neg.txt', path, generate_name(c, m, n)), vpa_digits).';
      everything.coefficients.N = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_N.txt', path, generate_name(c, m, n)), vpa_digits);
      everything.coefficients.log_abs_N = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_log_abs_N.txt', path, generate_name(c, m, n)), vpa_digits);
      everything.coefficients.F = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_F.txt', path, generate_name(c, m, n)), vpa_digits);
      everything.coefficients.log_abs_F = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_log_abs_F.txt', path, generate_name(c, m, n)), vpa_digits);
      everything.coefficients.k1 = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_k1.txt', path, generate_name(c, m, n)), vpa_digits);
      everything.coefficients.log_abs_k1 = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_log_abs_k1.txt', path, generate_name(c, m, n)), vpa_digits);
      everything.coefficients.k2 = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_k2.txt', path, generate_name(c, m, n)), vpa_digits);
      everything.coefficients.log_abs_k2 = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_log_abs_k2.txt', path, generate_name(c, m, n)), vpa_digits);
      everything.coefficients.c2k = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_c2k.txt', path, generate_name(c, m, n)), vpa_digits).';
      everything.coefficients.log_abs_c2k = pro_read_numerical_to_vpa(sprintf('%s/pro_%s_log_abs_c2k.txt', path, generate_name(c, m, n)), vpa_digits).';
    else
      everything.coefficients.lambda = load(sprintf('%s/pro_%s_lambda.txt', path, generate_name(c, m, n)));
      everything.coefficients.log_abs_lambda = load(sprintf('%s/pro_%s_log_abs_lambda.txt', path, generate_name(c, m, n)));
      everything.coefficients.dr = load(sprintf('%s/pro_%s_dr.txt', path, generate_name(c, m, n))).';
      everything.coefficients.log_abs_dr = load(sprintf('%s/pro_%s_log_abs_dr.txt', path, generate_name(c, m, n))).';
      everything.coefficients.dr_neg = load(sprintf('%s/pro_%s_dr_neg.txt', path, generate_name(c, m, n))).';
      everything.coefficients.log_abs_dr_neg = load(sprintf('%s/pro_%s_log_abs_dr_neg.txt', path, generate_name(c, m, n))).';
      everything.coefficients.N = load(sprintf('%s/pro_%s_N.txt', path, generate_name(c, m, n)));
      everything.coefficients.log_abs_N = load(sprintf('%s/pro_%s_log_abs_N.txt', path, generate_name(c, m, n)));
      everything.coefficients.F = load(sprintf('%s/pro_%s_F.txt', path, generate_name(c, m, n)));
      everything.coefficients.log_abs_F = load(sprintf('%s/pro_%s_log_abs_F.txt', path, generate_name(c, m, n)));
      everything.coefficients.k1 = load(sprintf('%s/pro_%s_k1.txt', path, generate_name(c, m, n)));
      everything.coefficients.log_abs_k1 = load(sprintf('%s/pro_%s_log_abs_k1.txt', path, generate_name(c, m, n)));
      everything.coefficients.k2 = load(sprintf('%s/pro_%s_k2.txt', path, generate_name(c, m, n)));
      everything.coefficients.log_abs_k2 = load(sprintf('%s/pro_%s_log_abs_k2.txt', path, generate_name(c, m, n)));
      everything.coefficients.c2k = load(sprintf('%s/pro_%s_c2k.txt', path, generate_name(c, m, n))).';
      everything.coefficients.log_abs_c2k = load(sprintf('%s/pro_%s_log_abs_c2k.txt', path, generate_name(c, m, n))).';
    endif
	catch
		fprintf('can''t open one or more of the coefficients...\n');
	end_try_catch
	try
    % reading PSWF function
    % NOTE : 'pro_sphwv' outputs in double precision - no need for 'vpa'
		everything.S1 = struct();
		unzip(sprintf('%s/pro_%s_%s.zip', path, generate_name(c, m, n), name), path);

		temp = load(sprintf('%s/pro_%s_S1.txt', path, generate_name(c, m, n))).';
		everything.S1.eta = temp(2, :);
		[ ...
		everything.S1.eta, ix ...
		] = sort(everything.S1.eta);
		everything.S1.S1_1 = temp(3, ix);
		everything.S1.S1p_1 = temp(4, ix);
		everything.S1.S1_2 = temp(5, ix);
		everything.S1.S1p_2 = temp(6, ix);
		everything.S1.S1_log_abs_difference = temp(7, ix);
		everything.S1.S1p_log_abs_difference = temp(8, ix);
		everything.S1.S1 = everything.S1.S1_1;
		everything.S1.S1p = everything.S1.S1p_1;
	catch
		fprintf('can''t open S1...\n');
  end_try_catch

  % compute via stable ratios "Prolate spheroidal wavefunctions, quadrature and interpolation"
  %       H Xiao, V Rokhlin and N Yarvin :
  %   for n = 0 use 1/fc(0)*int^{1}_{-1} fc(x)dx with a Simpson's rule
  %   for n > 1 use stable ratio scheme from:

  pswf = everything.S1.S1;
  pswf_grid = linspace(-1,1,length(pswf));
  if (n == 0)
    pswf_at_zero = interp1(pswf_grid, pswf, 0.0, "linear");
    pswf_int = IntegrateSimpson1D(pswf, pswf_grid(2)-pswf_grid(1));
    fc_eigval = pswf_int/pswf_at_zero;
  else
    % get previous PSWF and its eigenvalue
    pswf_struct_prev = pro_open_everything(saved, c, 0, n-1);
    pswf_prev = pswf_struct_prev.S1.S1;
    fc_eigval_abs_prev = abs(pswf_struct_prev.coefficients.mu);

    % get derivatives of pswfs
    d_pswf = diff(pswf)/(pswf_grid(2)-pswf_grid(1));
    d_pswf = interp1(pswf_grid(1:end-1), d_pswf, pswf_grid, "linear", "extrap");
    d_pswf_prev = diff(pswf_prev)/(pswf_grid(2)-pswf_grid(1));
    d_pswf_prev = interp1(pswf_grid(1:end-1), d_pswf_prev, pswf_grid, "linear", "extrap");

    % get nominator and denominator ratio
    pswf_int_nom_abs = abs(IntegrateSimpson1D(d_pswf_prev.*pswf, pswf_grid(2)-pswf_grid(1)));
    pswf_int_den_abs = abs(IntegrateSimpson1D(d_pswf.*pswf_prev, pswf_grid(2)-pswf_grid(1)));

    % compute new eigenvalue
    fc_eigval_abs = fc_eigval_abs_prev*sqrt(pswf_int_nom_abs/pswf_int_den_abs);
    fc_eigval = (1i^n)*fc_eigval_abs;
  endif
  everything.coefficients.mu = fc_eigval;
  everything.coefficients.mu

	delete(sprintf('%s/pro_%s_*.txt', path, generate_name(c, m, n)));
end
