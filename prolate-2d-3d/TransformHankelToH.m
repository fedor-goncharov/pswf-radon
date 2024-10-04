% transform Hankel transform to H-data
% see formulas (2.5), (2.7) in the original article

function [h_data, h_grid] = TransformHankelToH(hankel_data, bdwdth_grid, order, integer_order_flag=true)
  % params:
  %   hankel_data : Hankel integral transform array
  %   bdwdth_grid : parameter grid of transforms
  % Note:
  %   no global interpolation
  %   it assumed that bdwdth_grid(1) = 0.0
  %   to avoid zero divison at zero coordinate we use spline interpolation
  % Return:
  %   h_data : H-transform of original data
  %   h_grid : uniform grid in [-1,1]

  % set symmetric grid (note 0.0 at center) and symmetric data
  h_grid = [(-1.)*fliplr(bdwdth_grid(2:end)), bdwdth_grid];
  h_grid /= bdwdth_grid(end); % normalized grid on [-1,1]

  % symmetrize, scale and normalize data
  scale = [fliplr(bdwdth_grid(2:end)), 1.0, bdwdth_grid(2:end)];
  if integer_order_flag
    scale = 1./sqrt(scale);
    n = order;
  else
    scale = 1./scale;
    n = order-0.5;
  endif
  scale(1:length(bdwdth_grid)-1) *= (-1)^n;
  h_data = [fliplr(hankel_data(2:end)), hankel_data] .* scale;

  % interpolate at zero to avoid zero division
  no_zero_mask = true(1, length(h_grid));
  no_zero_mask(length(bdwdth_grid)) = false;
  h_data(length(bdwdth_grid)) = interp1(h_grid(no_zero_mask), h_data(no_zero_mask), 0.0, "spline");
endfunction
