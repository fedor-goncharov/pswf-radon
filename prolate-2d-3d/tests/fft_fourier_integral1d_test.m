% test on Fourier integral in 1d

test_func = @(x) (abs(x) < 0.5);

ngrid = 512;
grid = linspace(-1,1,ngrid + 1)(1:end-1);
Fs = 1./(grid(2) - grid(1)); % sampling rate (Nyquist frequency)
Fs_grid = linspace(-Fs/2, Fs/2, ngrid + 1)(1:end-1); % endpoint = false
test_func_sampled = arrayfun(test_func, grid);


% analytic Fourier transform
ft_func_real = (@(t) quad(@(x) test_func(x)*cos(-2*pi*x*t), -1, 1));
ft_func_imag = (@(t) quad(@(x) test_func(x)*sin(-2*pi*x*t), -1, 1));
ft_func = @(t) ft_func_real(t) + 1i*ft_func_imag(t);
ft_analytic = arrayfun(ft_func, Fs_grid);


% DFT approach
ft_dft =(1./Fs)*exp(-2*pi*1i*Fs_grid).*fftshift(fft(test_func_sampled));

figure, plot(Fs_grid, real(ft_analytic));
hold on;
plot(Fs_grid, real(ft_dft));

norm(ft_dft - ft_analytic)


% CONCLUSION: DFT approximates analytic Fourier integral
%             endpoint is CRUCIAL for precision
