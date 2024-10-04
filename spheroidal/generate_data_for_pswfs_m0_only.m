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

function generate_data_for_pswfs_m0_only(c, nmax=29, nintervals=1024)
  m = 0
  for n = m : m + nmax
		pro_calculate_everything(sprintf('%s/sphwv', pwd()), 2000, 500, c, m, n, 10, -300, 10, -300, 10, -300);
		everything = pro_load_everything('data', c, m, n, '');
		pro_save_everything('saved', everything);
	end

	for n = m : m + nmax
		n % trace current n

		pro_calculate_functions('', sprintf('%s/sphwv', pwd()), 2000, 500, c, m, n, 2.0 / nintervals, 'eta', -7, 9.0, 1.0 / 64.0, 'x', 17);
		everything = pro_load_everything_with_eigvalues('data', 'saved', c, m, n, '');
		pro_save_everything('saved', everything);
	end
end
