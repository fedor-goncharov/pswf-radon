% handy function to 'imshow' directly from radial data


function plot2dfrom1d(f_data, f_grid, title_str="", symmetrize_grid = true)

  grid = f_grid;
  if symmetrize_grid == true
    grid = [-fliplr(f_grid)(1:end-1) f_grid];
  endif
  [YY, XX] = ndgrid(grid, grid);
  RR = sqrt(YY.^2 + XX.^2);

  vals2d = interp1(f_grid, f_data, RR, "linear", 0.0);
  figure, imagesc(vals2d, [min(min(vals2d)), max(max(vals2d))]), colorbar, title(title_str);
endfunction
