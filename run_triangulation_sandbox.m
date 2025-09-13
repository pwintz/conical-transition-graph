clear;
dim = 3;
% dim = 2;
switch dim
  case 2
    conical_partition = ConicalPartition.fromNumSlices(16);
  case 3
    sphere_grid = sphereGrid3d(nLinesOfLatitude=10, nLinesOfLongitude=10);
    conical_partition = ConicalPartition(sphere_grid);
end

tri = conical_partition.triangulation;


% ╭─────────────────────────────────────────────────╮
% │             Get the number of cones             │
% ╰─────────────────────────────────────────────────╯
n_cones_x_verts_per_cone = tri.size;
n_cones        = n_cones_x_verts_per_cone(1);
verts_per_cone = n_cones_x_verts_per_cone(2);

pwintz.assertions.assertEqual(n_cones, size(tri.ConnectivityList, 1));

% ╭────────────────────────────────────────────────────────╮
% │             Get the cone contaning a point             │
% ╰────────────────────────────────────────────────────────╯
point = 0.1*ones(dim, 1);% Must be inside the convex hull of the vertices.
cone_ndx = tri.pointLocation(point');

% ╭────────────────────────────────────╮
% │             Get Origin             │
% ╰────────────────────────────────────╯
origin = zeros(1, dim);
origin_ndx = tri.nearestNeighbor(origin);

% ╭──────────────────────────────────────────────╮
% │             Get Vertices of Cone             │
% ╰──────────────────────────────────────────────╯
vertex_ndxs = tri.ConnectivityList(cone_ndx, :);
vertices    = tri.Points(vertex_ndxs, :)';

% ╭───────────────────────────────────────────────────────────╮
% │             Get cones adjacent to each vertex             │
% ╰───────────────────────────────────────────────────────────╯
adjacent_cone_ndxs = tri.vertexAttachments();

% ╭────────────────────────────────────────────╮
% │             Get adjacent cones             │
% ╰────────────────────────────────────────────╯
% "neighbors" function returns 'nan' if there is an external face.
nb_ndxs_with_nan = tri.neighbors();
nb_ndxs = nb_ndxs_with_nan(~isnan(nb_ndxs_with_nan));

% ╭────────────────────────────────────────────────────────╮
% │             Get vertices connected by edge             │
% ╰────────────────────────────────────────────────────────╯
vertex_ndx = vertex_ndxs(1);
[edges_of_vertex, ~] = find(vertex_ndx == tri.edges());

pwintz.strings.format("Edge indices connected to vertex %d:\n\t%d", vertex_ndx, edges_of_vertex)

% ╭──────────────────────────────────────────╮
% │             Get Cone Centers             │
% ╰──────────────────────────────────────────╯
centers = tri.incenter();

% ╭──────────────────────────────────────────╮
% │             Get Cone Normals             │
% ╰──────────────────────────────────────────╯
cone_hrep = HalfspaceRepresentation.fromConicalHull(vertices);
% [A_ineq, b_ineq, A_eq, b_eq] = polyhedron.vert2lcon(vertices(:, 2:end)')
% h.A_ineq



% [facets_vertex_ndxs, boundary_points] = tri.freeBoundary()
% % facets_vertex_ndxs = tri.freeBoundary();
% pwintz.strings.format("Facets (%z): %s, ", facets_vertex_ndxs, facets_vertex_ndxs)

% boundary_triangulation = triangulation(facets_vertex_ndxs, boundary_points)
% boundary_triangulation.edges()

% pwintz.plots.plotPoints(boundary_points')


pwintz.plots.namedFigure("Conical Partition Cones");
clf();
xlim(1.2*[-1, 1]);
ylim(1.2*[-1, 1]);
zlim(1.2*[-1, 1]);
view(3); % set the default 3D view, AZ = 0, EL = 90.
axis square;
hold on;
legend();
% conical_partition.plot()


disp("Plotting cones.");
for cone_ndx = conical_partition.cone_indices
  conical_partition.getCone(cone_ndx);
  conical_partition.plotCone(cone_ndx);
  % drawnow();
  % pause(1);
end 

disp("Plotting cone boundaries.");
for bnd_ndx = conical_partition.boundary_indices
  % conical_partition.getCone(cone_ndx)
  conical_partition.plotBoundary(bnd_ndx);
  % drawnow();
  % pause(0.01);
end 

disp("Plotting cone sphere approximations.");
for cone_ndx = conical_partition.cone_indices
  % conical_partition.getCone(cone_ndx)
  conical_partition.plotOverapproximationOfUnitSphereInCone(cone_ndx);
  % drawnow();
  % pause(0.1);
end 

return

disp("Plotting cone boundaries sphere approximations.");
for bnd_ndx = conical_partition.boundary_indices
  % conical_partition.getCone(cone_ndx)
  conical_partition.plotOverapproximationOfUnitSphereInBoundary(bnd_ndx);
  drawnow();
  % pause(1);
end 



return

if dim == 2
  pwintz.plots.namedFigure("2D Triangularization");
  clf();
  triplot(tri, "displayName", "Triangularization");
  hold on;
  pwintz.plots.plotPoints(centers', plotArgs={"DisplayName", "Centers"});
  legend();
end

if dim == 3

  pwintz.plots.namedFigure("Facets of 3D Triangularization");
  clf();
  clf();
  xlim(1.3*[-1, 1]);
  ylim(1.3*[-1, 1]);
  zlim(1.3*[-1, 1]);
  axis square;
  hold on;
  trisurf(facets, points(:,1),points(:,2),points(:,3), 'FaceColor','cyan','FaceAlpha',0.2);
  pwintz.plots.plotPoints(centers', plotArgs={"DisplayName", "Centers"});

  % ╭───────────────────────────────────────────╮
  % │             Get Boundary Cone             │
  % ╰───────────────────────────────────────────╯

  pwintz.plots.namedFigure("Boundary cones of a cone in 3D Triangularization");
  clf();
  xlim(1.3*[-1, 1]);
  ylim(1.3*[-1, 1]);
  zlim(1.3*[-1, 1]);
  axis square;
  pwintz.plots.plotPoints(sphere_grid, plotArgs={"Color", "black", "Marker", ".", "DisplayName", "Sphere"});
  hold on;
  pwintz.plots.plotPoints([0; 0; 0], plotArgs={"Color", "black", "Marker", "o", "DisplayName", "Origin"});
  pwintz.plots.plotPoints(vertices, plotArgs={"Color", "red", "Marker", "*", "DisplayName", pwintz.strings.format("vertex\\_ndxs %d", vertex_ndxs)});
  for nb_cone_ndx = nb_ndxs
    % cone_ndx
    % nb_cone_ndx
    % vertices
    nb_vertex_ndxs = tri.ConnectivityList(nb_cone_ndx, :);
    nb_vertex_ndxs_wo_origin = setdiff(nb_vertex_ndxs, origin_ndx);
    nb_vertices    = tri.Points([nb_vertex_ndxs_wo_origin, nb_vertex_ndxs_wo_origin(1)], :)';

    bnd_vertex_ndxs = intersect(nb_vertex_ndxs, vertex_ndxs);
    bnd_vertices    = tri.Points([bnd_vertex_ndxs, bnd_vertex_ndxs(1)], :)';
    

    pwintz.plots.plotPoints(nb_vertices, plotArgs={"Marker", "s", "LineStyle", "-", "DisplayName", pwintz.strings.format("nb\\_vertices %d", nb_vertex_ndxs)});
    pwintz.plots.plotPoints(bnd_vertices, plotArgs={"LineStyle", ":", "LineWidth", 6, "DisplayName", "bnd\_vertices"});
    

    patch("XData", bnd_vertices(1, :), "YData", bnd_vertices(2, :), "ZData", bnd_vertices(3, :), "FaceAlpha", 0.1, HandleVisibility='off');
    legend();

    bnd_cone_hrep = HalfspaceRepresentation.fromConicalHull(bnd_vertices);

  end


end


% ╭───────────────────────────────────────────────────────╮
% │  ╭─────────────────────────────────────────────────╮  │
% │  │             Plotting Triangulations             │  │
% │  ╰─────────────────────────────────────────────────╯  │
% ╰───────────────────────────────────────────────────────╯
% ╭────────────────────────────────────────╮
% │             Plotting in 3D             │
% ╰────────────────────────────────────────╯
% To plot a triangulation in 3D as a collection of tetrahedrons:
tetramesh(tirang);