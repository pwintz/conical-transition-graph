
pwintz.plots.namedFigure("Eigenvalue Test");
clf;
xlim(10*[-1, 1]);
ylim(10*[-1, 1]);
axis square;
hold on;

n = 200;
n_real = 0;
n_complex = 0;
for i = 1:100
  A = randn(n);
  eigenvals = eig(A);
  
  n_real = n_real + sum(imag(eigenvals) == 0);
  n_complex = n_complex  + sum(imag(eigenvals) ~= 0)/2;
  plot(real(eigenvals), imag(eigenvals), ".");
end

fprintf('   Percent real: %.3g%%\n', 100 * n_real / (n_real + n_complex));
fprintf('Percent complex: %.3g%%\n', 100 * n_complex / (n_real + n_complex));