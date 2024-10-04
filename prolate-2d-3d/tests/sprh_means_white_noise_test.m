
domain_grid_size = 256  ;
lin = linspace(-1,1, domain_grid_size + 1)(1:end-1);
[XX, YY] = meshgrid(lin, lin);
dx = lin(2) - lin(1);
cell_volume = dx*dx;

phi_grid_size = 256;
phi_grid = linspace(0, 2*pi, phi_grid_size + 1)(1:end-1);
dphi = phi_grid(2) - phi_grid(1);

xs = vec(cos(phi_grid));
ys = vec(sin(phi_grid));

r_grid_size = 64;
r_grid = linspace(0, 1, r_grid_size);
dr = r_grid(2) - r_grid(1);

nsamples = 2000;
cont_sph_means = zeros(r_grid_size, nsamples);

sigma = 1.0

for n=1:nsamples
  n
  gfield = sigma*randn(domain_grid_size, domain_grid_size)/sqrt(cell_volume); % discrete appr. to white noise
  sph_means = arrayfun(@(r) trapz(interp2(XX, YY, gfield, r*xs, r*ys))*dphi, vec(r_grid)).*sqrt(vec(r_grid));
  cont_sph_means(:, n) = sph_means;
endfor

xs_xy_slice = vec(cos(phi_grid));
ys_xy_slice = vec(sin(phi_grid));

xs_xyz = repmat(xs_xy_slice, [1 nsamples]);
ys_xyz = repmat(ys_xy_slice, [1 nsamples]);
zs_xyz = repmat([1:nsamples], [phi_grid_size 1]);
gfield = sigma/sqrt(cell_volume)*randn(domain_grid_size, domain_grid_size, nsamples);

[XX, YY, ZZ] = meshgrid(lin, lin, 1:nsamples);

cont_sph_means = zeros(r_grid_size, nsamples);
for ir=1:r_grid_size
  ir
  r = r_grid(ir);
  sph_means = trapz(interp3(XX, YY, ZZ, gfield, r*xs_xyz, r*ys_xyz, zs_xyz), 1)*dphi;
  cont_sph_means(ir, :) = sph_means(:);
endfor



%sph_means = arrayfun(@(r)

sqrt_rho_sph_means = cont_sph_means.*sqrt(vec(r_grid));
figure, imagesc(cov(sqrt_rho_sph_means')), colorbar;
colorbar;
figure, imagesc(cov(sqrt_rho_sph_means')), colorbar;

plot(1:256, sqrt_rho_sph_means(:, 128));
