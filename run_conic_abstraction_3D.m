% ╭────────────────────────────────────────────────────────────╮
% │  ╭──────────────────────────────────────────────────────╮  │
% │  │             3D Conic Abstraction Example             │  │
% │  ╰──────────────────────────────────────────────────────╯  │
% ╰────────────────────────────────────────────────────────────╯
% rng(1);
randomize_matrices = false;
randomize_matrices = true;
randomize_matrices = ~exist("Ac", "var") || ~exist("Ad", "var") || randomize_matrices;
if randomize_matrices
  Ac = randn(3);
  Ad = randn(3);
end
pwintz.assertions.assertSize(Ac, [3, 3]);
pwintz.assertions.assertSize(Ad, [3, 3]);

[Ac_eigvecs, Ac_eigvals] = eig(Ac, "vector");
if any(real(Ac_eigvals) >= 0) 
% for i = 1:size(Ac, 1);
  % real(Ac_eigvals(i))
  warning("Some of the eigenvalues of Ac, %s, are in the right-half of the complex plane.", mat2str(Ac_eigvals, 3));
end

flow_set_angles = [0, pi];
jump_set_angles = [pi-0.8, pi+0.0];
conic_abstraction = ConicAbstraction.fromUVSphere3D(...
nLinesOfLatitude =10, ...
nLinesOfLongitude=10);
% conic_abstraction = ConicAbstraction.fromUVSphere3D(...
%   "flowMapMatrix", Ac, ...
%   "jumpMapMatrix", Ad ...
% );


% ╭────────────────────────────────────╮
% │  ╭──────────────────────────────╮  │
% │  │             Plot             │  │
% │  ╰──────────────────────────────╯  │
% ╰────────────────────────────────────╯
% ╭────────────────────────────────────────────────────╮
% │             Plot the conic abstraction             │
% ╰────────────────────────────────────────────────────╯
pwintz.plots.namedFigure("3D Conic Abstraction");
clf();
xlim(1.2*[-1, 1]);
ylim(1.2*[-1, 1]);
zlim(1.2*[-1, 1]);
view(3);
axis square;
hold on;
conic_abstraction.plotCones();