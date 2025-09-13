% rng(1);
randomize_matrices = false;
randomize_matrices = true;
randomize_matrices = ~exist("Ac", "var") || ~exist("Ad", "var") || ~all(size(Ac) == [2, 2]) || randomize_matrices;
if randomize_matrices
  Ac = randn(2, 2);
  Ad = randn(2, 2);
end

Ac = [-1, -0.6; 
       5, -1]
Ad = [
  -1, 1;
  -0.5, -5;
]

% ! These matrices result in a nice alignment of C, D, and G(D). 
% Ac = [
% 	  -0.1511,   -0.3642;
% 	   0.8928,   -0.8262];
% Ad = [-1.5180,   -3.0722;
%       -1.0749,    0.5214];

% ! This combination of matrices looks like it would produce stable flows but actually is unstable. 
% Ac =
%    -0.2730   -0.4809
%     1.5763    0.3275
% Ad =
%     0.0799    0.4115
%    -0.9485    0.6770

[Ac_eigvecs, Ac_eigvals] = eig(Ac, "vector");
if any(real(Ac_eigvals) >= 0) 
% for i = 1:size(Ac, 1);
  % real(Ac_eigvals(i))
  warning("Some of the eigenvalues of Ac, %s, are in the right-half of the complex plane.", mat2str(Ac_eigvals, 3));
end

% conical_partition = ConicalPartition.fromNumSlices(10);
flow_set_angles = [0, pi];
jump_set_angles = [pi-0.2, pi+0.0];
conic_abstraction = ConicAbstraction.fromAngles2D(...
  "flowMapMatrix", Ac, ...
  "jumpMapMatrix", Ad, ...
  "flowSetAngles", flow_set_angles, ...
  "jumpSetAngles", jump_set_angles, ...
  "maxStateConeAngle", 2*pi/4, ...
  "maxDerivativeConeAngle", 2*pi/100 ...
);
conical_partition = conic_abstraction.conical_partition;
flow_graph        = conic_abstraction.flow_transition_graph;
jump_graph        = conic_abstraction.jump_transition_graph;
ctg               = conic_abstraction.ctg;
tri               = conical_partition.triangulation;

pwintz.plots.namedFigure("flow transition graph");
clf();
xlim(1.3*[-1, 1]);
ylim(1.3*[-1, 1]);
axis square;
hold on;
conic_abstraction.plotCones();
conic_abstraction.flow_transition_graph.plot();
conic_abstraction.jump_transition_graph.plot();

conic_abstraction.hasStableFlows()


% ╭───────────────────────────────────────────────────────────────────────╮
% │             Test why one edge was missing from flow graph             │
% ╰───────────────────────────────────────────────────────────────────────╯
% vertex_14 = conical_partition.getVertex(14);
% cone_14 = conical_partition.getCone(14);
% [reachable_set_in_cone, reachable_set] = flow_graph.computeReachableSet(14, vertex_14);

% vertex_13 = conical_partition.getVertex(13)
% cone_13 = conical_partition.getCone(13)
% [reachable_set_in_cone, reachable_set] = flow_graph.computeReachableSet(13, vertex_13)
% 
% figure(2); clf
% clf();
% xlim(1.2*[-1, -0.8]);
% ylim(0.2*[-1, 1]);
% axis square;
% hold on;
% plot(1e4*cone_13, "FaceColor", "red", "FaceAlpha", 0.5)
% plot(reachable_set_in_cone, "FaceColor", "green", "FaceAlpha", 0.1);
% R_from_13 = conic_abstraction.flow_transition_graph.directly_reachable_sets_from_vertex_through_cone(14, :);
% R = R_from_13{13};
% 
% conic_abstraction.flow_transition_graph.getEdgesFromVerticesToVertices();
% this_should_be_an_edge = conic_abstraction.flow_transition_graph.getEdgesFromVertexToVertex(14, 15);
% assert(~isempty(this_should_be_an_edge), 'We expect there to be an edge from vertex 14 to vertex 15');
% this_should_be_an_edge = conic_abstraction.flow_transition_graph.getEdgesFromVertexToVertex(13, 14);
% assert(~isempty(this_should_be_an_edge), 'We expect there to be an edge from vertex 13 to vertex 14');

conic_abstraction.flow_transition_graph.getEdgesFromVerticesToVertices();
conic_abstraction.flow_transition_graph.getEdgesFromVerticesToCones();
conic_abstraction.flow_transition_graph.getEdgesFromConesToVertices();
conic_abstraction.flow_transition_graph.getEdgesFromConesToCones();
% conic_abstraction.flow_transition_graph.getVertexSubgraph().Nodes
% conic_abstraction.flow_transition_graph.getVertexSubgraph().Edges
% conic_abstraction.flow_transition_graph.getConeSubgraph().Nodes
% conic_abstraction.flow_transition_graph.getConeSubgraph().Edges
% vertex_subgraph = conic_abstraction.flow_transition_graph.getVertexSubgraph()
% contracted_flow_transition_graph = conic_abstraction.flow_transition_graph.contractEdgesToBetweenCones(conic_abstraction.jump_set_image_cone_ndxs, conic_abstraction.jump_set_cone_ndxs);
% disp(contracted_flow_transition_graph);

% ctg = TransitionGainDigraph.union(jump_graph, contracted_flow_transition_graph);
% [cycle_min_gains, cycle_max_gains, cycles_nodes, cycles_edges] =  ctg.getCycleGains();

% for start_cone_ndx = conic_abstraction.jump_set_image_cone_ndxs
%   for start_vertex_ndx = find(conic_abstraction.can_flow_from_cone_to_vertex(start_cone_ndx, :))
%     for end_cone_ndx = conic_abstraction.jump_set_cone_ndxs
%       for end_vertex_ndx = find(conic_abstraction.can_flow_from_vertex_into_cone(end_cone_ndx, :))
%       
%         % fprintf('Checking if there is a path from vertex %d to vertex %d\n', start_vertex_ndx, end_vertex_ndx);
%         [paths, path_edge_ndxs] = conic_abstraction.flow_transition_graph.getPathsBetweenVertices(start_vertex_ndx, end_vertex_ndx);
%         if ~isempty(paths)
%           % paths{1}
%           % path_edge_ndxs{1}
%           if length(paths) > 1
%             fprintf('There are multiple (%d) paths from cone %d to vertex %d to vertex %d to cone %d. They are.\n', length(paths), start_cone_ndx, start_vertex_ndx, end_cone_ndx, end_vertex_ndx);
%             disp(paths);
%           else
%             paths
%             path_edge_ndxs
%             edges  = conic_abstraction.flow_transition_graph.getEdgeRow([path_edge_ndxs{:}])
%             edge_weights = edges.MaxGain'
%             path_weight = prod(edge_weights)
%             fprintf('There is one path from cone %2d to vertex %2d to vertex %2d to cone %2d.\t', start_cone_ndx, start_vertex_ndx, end_cone_ndx, end_vertex_ndx);
%             fprintf('The path length was %2d edges and had a weight of %s\n', length(path_edge_ndxs{1}), ctg.utils.num2strNear1(path_weight));
%           end
%         end
%       end
%     end
%   end
% end

pwintz.plots.namedFigure("Reachable Cones from vertices");
clf();
xlim(1.2*[-1, 1]);
ylim(1.2*[-1, 1]);
axis square;
hold on;
% conic_abstraction.plotVertices();

% flow_graph.Edges
% matgeom_flow_graph = struct( ...
%   "nodes", conical_partition.vertices', ...
%   "edges", flow_graph.Edges.VertexNdxs ...
% )
% drawGraph(matgeom_flow_graph)

pwintz.plots.plotUnitCircle();
pwintz.plots.plotLinearVectorField(Ac);
conic_abstraction.plotCones();


conic_abstraction.ctg.plot();
conic_abstraction.plotConesReachableFromVertices();
% conic_abstraction.plotVerticesReachableFromCones();
% conic_abstraction.plotFlowGraph();
% conic_abstraction.plotJumpGraph();

% contracted_flow_transition_graph.plot();

pwintz.strings.format("      n_flow_set_cones: %d\t       flow_set_cone_ndxs: %d", conic_abstraction.n_flow_set_cones, conic_abstraction.flow_set_cone_ndxs)
pwintz.strings.format("      n_jump_set_cones: %d\t       jump_set_cone_ndxs: %d", conic_abstraction.n_jump_set_cones, conic_abstraction.jump_set_cone_ndxs)
pwintz.strings.format("n_jump_set_image_cones: %d\t jump_set_image_cone_ndxs: %d", conic_abstraction.n_jump_set_image_cones, conic_abstraction.jump_set_image_cone_ndxs)

if  conic_abstraction.hasStableFlows()
  disp("The flow component of the system is stable.");
else
  disp("The flow component of the system may be unstable.");
end % End of if "conic_abstraction.hasStableFlows()"

if conic_abstraction.is_origin_asymptotically_stable
  fprintf("The origin of the system is stable!\n");
else 
  fprintf("The system may be unstable!\n");
end
return % TODO: Remove


for cone_ndx = conical_partition.cone_indices
  circle_overapproximation = conical_partition.getUnitSphereOverApproximationInCone(cone_ndx);
  plot(circle_overapproximation, "FaceColor", "none", "EdgeColor", "black");
end


% for reachable_set = [conic_abstraction.reachable_sets_from_unit_sphere_in_cone{:}]
%   plot(reachable_set, "FaceColor", [0, 0.6, 0.7], "EdgeColor", "none", "FaceAlpha", 0.1);
% end % End of for block


for restricted_reachable_set = [conic_abstraction.restricted_reachable_sets_from_unit_sphere_in_cone{:}]
  plot(restricted_reachable_set, "FaceColor", [0, 0.6, 0.7], "EdgeColor", "none", "FaceAlpha", 0.6);
end % End of for block

for vertex_ndx = conical_partition.vertex_indices
  for reachable_set_from_vertex = [conic_abstraction.directly_reachable_sets_from_vertices{vertex_ndx, :}]
    plot(reachable_set_from_vertex, "FaceColor", [0.4, 0.0, 0.0], "EdgeColor", "none", "FaceAlpha", 0.9);
  end
end

return % TODO: Remove

angles = pwintz.arrays.range(0, 2*pi, ...
   "n_values", 20, "includeEnd", false);
conical_partition = ConicalPartition.fromAngles2D(angles);

graph = conical_partition.mesh_graph;
vertex_table = conical_partition.vertex_table;
% cone_table   = conical_partition.cone_table;
vertex_table.('position');
cell2mat(vertex_table{1, "adjacent_cones_ndxs"});
cell2mat(vertex_table{2, "adjacent_cones_ndxs"});

vertex_table.adjacent_cones_ndxs{1};

face_vertices = [graph.nodes(graph.faces(:, 1), :), graph.nodes(graph.faces(:, 2), :), graph.nodes(graph.faces(:, 3), :)];

[polygon] = grFaceToPolygon(graph, 1);

% ╭────────────────────────────────────────────────────────────╮
% │             Draw Condes Adjacent to Boundaries             │
% ╰────────────────────────────────────────────────────────────╯
for ray_ndx = conical_partition.ray_indices
  ray = conical_partition.getRay(ray_ndx);
  adjacent_cone_ndxs = conical_partition.getConesAdjacentToRay(ray_ndx);
  figure(1);
  clf();
  xlim(1.0*[-1, 1]);
  ylim(1.0*[-1, 1]);
  axis square;
  hold on;
  pwintz.plots.plotVector2(ray, plotArgs={"LineWidth", 12})
  adjacent_cones = conical_partition.getCone(adjacent_cone_ndxs)
  testCase.assertNumElements(adjacent_cone_ndxs, 2, "Each rays whould have 2 adjacent_cone_ndxs");
  testCase.assertNumElements(adjacent_cones, 2, "Each rays whould have 2 adjacent_cones");
  for cone = adjacent_cones
    assert(cone.n_rays == 2)
    cone.plot();
    % ray.plot()
    drawnow();
    pause(1);
    % testCase.assertTrue(cone.contains(ray));
    
  end
end


% vertex_table
% cone_table
return % TODO: Remove
% % conical_partition
% % conical_partition.mesh_graph
% [nodes, edges] = gabrielGraph([conical_partition.vertices'; [0, 0]]);
% 
% % origin_ndx = find(nodes(:, 1) == 0 & nodes(:, 2) == 0);
% origin_ndx =  pwintz.arrays.findRowIn([0, 0], nodes)
% not_origin_ndxs = find(1:size(nodes, 1) ~= origin_ndx)
% 
% % ⋘──────── Construct faces ────────⋙
% 
% faces_ndx_col_1 = not_origin_ndxs'
% faces_ndx_col_2 = circshift(not_origin_ndxs,1)'
% faces_ndx_col_3 = origin_ndx * ones(size(not_origin_ndxs'))
% faces = [col_1, col_2, col_3]
% 
% 
figure(1);
clf;
% drawGraph(nodes, edges);
axis square;
xlim(1.2*[-1, 1]);
ylim(1.2*[-1, 1]);
hold on;
% 
% grAdjacentNodes(edges, 2);
% 
% graph = struct( ...
%   "nodes", nodes, ...
%   "edges", edges, ...
%   "faces", faces ...
% )

colors = {'red', 'blue', 'green'};
for i=1:length(graph.faces)
  [polygon] = grFaceToPolygon(graph, i);
  fillPolygon(polygon, colors{mod(i, numel(colors)) + 1});
  hold on;
end

% graph.nodes
% graph.edges
% graph.faces

% fillGraphFaces(graph, )

% all_verts = [conical_partition.vertices'; [0, 0]];
% [NODES, EDGES, FACES] = voronoi2d([all_verts; 2*conical_partition.vertices'])
% 
% figure(2);
% clf();
% drawGraph(NODES, EDGES);
% hold on
% pwintz.plots.plotUnitCircle();