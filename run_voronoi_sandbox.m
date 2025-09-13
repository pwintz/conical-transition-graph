clear
% ╭───────────────────────────────────────────╮
% │             Generate Vertices             │
% ╰───────────────────────────────────────────╯

dimension = 2;

switch dimension
  case 2
    thetas  = linspace(0, 2*pi, 10);
    
    polar2cartesian = @(r, theta) [
      r .* cos(theta);
      r .* sin(theta);
    ];
    polar2cartesian(1, 2)
    
    nodes = pwintz.arrays.mapRows(@(theta) polar2cartesian(1, theta)', thetas');
    pwintz.assertions.assertAlmostEqual(vecnorm(nodes'), 1);

  case 3
    % ⋘──────── Generate vertices from angles ────────⋙
    phis = linspace(0, 2*pi, 10);
    thetas  = linspace(0, pi, 10);

    [Phi, Theta] = meshgrid(phis, thetas);
    Phi = reshape(Phi, [], 1);
    Theta = reshape(Theta, [], 1);

    spherical2cartesian = @(r, theta, phi) [
      r .* sin(theta) .* cos(phi);
      r .* sin(theta) .* sin(phi);
      r .* cos(theta);
    ];


    nodes = pwintz.arrays.mapRows(@(theta, phi) spherical2cartesian(1, theta, phi)', Theta, Phi);
    pwintz.assertions.assertAlmostEqual(vecnorm(nodes'), 1);

  % ⋘──────── Generate vertices random  ly ────────⋙
  % nodes = 2*rand(9, 3) - 1
  % nodes = (nodes' ./ vecnorm(nodes'))'
  otherwise
    error("Unexpected case.");
end

% ⋘──────── Add origin ────────⋙
origin = 0*nodes(1, :);
nodes = [nodes; origin];
origin_ndxs = pwintz.arrays.findRowIn(origin, nodes);

% ⋘──────── Remove duplicates ────────⋙
nodes = pwintz.arrays.uniqueRows(nodes, tolerance=1e-8);

% ╭────────────────────────────────────────╮
% │             Generate edges             │
% ╰────────────────────────────────────────╯

% [nodes, edges]  = gabrielGraph(nodes);

triang = delaunayTriangulation(nodes)
triang.ConnectivityList

figure(1); 
clf();
xlim(1.1*[-1, 1]);
ylim(1.1*[-1, 1]);
zlim(1.1*[-1, 1]);
axis square;
hold on;

centers = triang.incenter();
for i = 1:size(triang.ConnectivityList, 1)
  ndxs = triang.ConnectivityList(i, :);
  points = nodes(setdiff(ndxs, []), :)
  points = [points; points(1, :)]
  
  center = centers(i, :);
  
  % plot3(points(:, 1), points(:, 2), points(:, 3), "-");
  % plot3(center(1), center(2), center(3), "*")
  plot(points(:, 1), points(:, 2), "-");
  plot(center(1), center(2), "*")
  drawnow
  pause(0.1)
end
return

figure(1);
clf

drawGraph(nodes, edges);

return
% [nodes, edges] = delaunayGraph(nodes);

% zero_node_ndx = pwintz.arrays.findRowIn([0, 0, 0], nodes)
% nodes = nodes(setdiff(1:end, zero_node_ndx), :)
% edges = edges(edges(:, 1) ~= zero_node_ndx, :)
% edges = edges(edges(:, 2) ~= zero_node_ndx, :)


% ╭────────────────────────────────────╮
% │             Find faces             │
% ╰────────────────────────────────────╯
% ⋘──────── Use cycles ────────⋙
% struct2table(struct(...
%   "do", value1, ...
%   "name2", value2  ...
% ));
edge_table = table(edges,'VariableNames',{'EndNodes'})

G = digraph(edge_table);

edge_table
G = digraph(edge_table);

figure(1);
clf;
% plot(G)


drawGraph(nodes, edges);


return






% ╭─────────────Voronoi─────────────╮
% │             Voronoi             │
% ╰─────────────Voronoi─────────────╯

[voronoi_nodes,C] = voronoin(nodes)
voronoi_nodes

voronoi_nodes = (voronoi_nodes' ./ vecnorm(voronoi_nodes'))'

pwintz.plots.namedFigure("Voronoi");
clf();
xlim(2*[-1, 1]);
ylim(2*[-1, 1]);
zlim(2*[-1, 1]);
axis square;
hold on;
plot3(nodes(:, 1), nodes(:, 2), nodes(:, 3), "*", "DisplayName", "original points");
hold on;
plot3(voronoi_nodes(:, 1), voronoi_nodes(:, 2), voronoi_nodes(:, 3), "*", "DisplayName", "vertices");


drawGraph(nodes, edges);

legend();