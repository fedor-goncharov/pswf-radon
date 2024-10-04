  % test harmonicY library for orthonormality


  phi_grid_size = 256;
  theta_grid_size = 256;
  phi_grid = linspace(0, 2*pi, phi_grid_size);
  theta_u_grid = linspace(-1, 1, theta_grid_size + 2)(2:end-1);
  theta_grid = acos(theta_u_grid);

  [PP, TT] = meshgrid(phi_grid, theta_grid); % 2D grid of angles (theta, phi)

  sph0 = harmonicY(0, 0, TT, PP, 'norm', true);
  sph2 = harmonicY(2, 0, TT, PP, 'norm', true);
  sph3 = harmonicY(3, 0, TT, PP, 'norm', true);

  sum(sph0(:).^2)*(2*pi)/(phi_grid_size-1)*(2/(theta_grid_size-1))
  sum(sph2(:).^2)*(2*pi)/(phi_grid_size-1)*(2/(theta_grid_size-1))
  sum(sph3(:).^2)*(2*pi)/(phi_grid_size-1)*(2/(theta_grid_size-1))
