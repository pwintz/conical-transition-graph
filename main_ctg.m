

% roll_array = linspace(0, 2*pi, 10)
% pitch_array = linspace(-pi, pi, 10)
% yaw_array = linspace(0, pi, 10);

% ╭─────────────────────────────────────╮
% │ ╭─────────────────────────────────╮ │
% │ │             Options             │ │
% │ ╰─────────────────────────────────╯ │
% ╰─────────────────────────────────────╯
clear();
% Create a random 2x2 matrix to define \dot x = A_c x.
A_cont = rand(2) - 0.5;
% Create a random 2x2 matrix to define x^+ = A_d x.
A_disc = rand(2) - 0.5;
n_slices = 40;


% ╭────────────────────────────────────────╮
% │ ╭────────────────────────────────────╮ │
% │ │             Given Data             │ │
% │ ╰────────────────────────────────────╯ │
% ╰────────────────────────────────────────╯

flow_set = ConvexPolyhedron.fromConvexHull([[-1; 0], [0; 0], [1; 0.5]]);
% assert(flow_set, "error message format string", error message parameters)

jump_set        = 1e3 * ConvexPolyhedron.fromConvexHull([[1; 0], [0; 0], [1; -0.01]]);
jump_set_image  = 1e3 * A_disc * jump_set;

% ⋘────────── Some interesting choices of A_cont ──────────⋙
% # This A_cont gives a very large slices in the state space. One of its eigenvalues is -0.007, so the problem likely arises from the nearly zero eigenvalue.
% A_cont =
%     0.4027   -0.0091
%     0.4448   -0.0107
% # This A_cont creates an interesting alignment, with the flow pointing nearly exactly the opposite direction as the state cone.
% A_cont = [
% 	  -0.3103,   -0.3524;
% 	  -0.0050,   -0.4450
% ]

% ╭────────────────────────────────────────────╮
% │ ╭────────────────────────────────────────╮ │
% │ │             CTG Parameters             │ │
% │ ╰────────────────────────────────────────╯ │
% ╰────────────────────────────────────────────╯

conical_partition = ConicalPartition(A_cont, "nSlices", n_slices);
[jump_set_nonempty_intersection_ndxs, jump_set_in_regions] = conical_partition.getStateRegionIntersections(jump_set);
[jump_set_image_nonempty_intersection_ndxs, jump_set_image_in_regions] = conical_partition.getStateRegionIntersections(jump_set_image);
jump_set_nonempty_intersection_ndxs
jump_set_image_nonempty_intersection_ndxs
% jump_set_image
assert(~isempty(jump_set_nonempty_intersection_ndxs),       "No regions were found to intersect with the jump set. This should not happen because the sets are cones, but is possible because we stored them as large--but bounded polytopes.");
assert(~isempty(jump_set_image_nonempty_intersection_ndxs), "No regions were found to intersect with the image of the jump set. This should not happen because the sets are cones, but is possible because we stored them as large--but bounded polytopes.");
return


conic_abstraction = ConicAbstraction(conical_partition, flow_set);
graph = conic_abstraction.graph;
% cyclebasis(graph)

return
% Get all of the fundamental cycles in the graph. 
[cycles, edgecycles] = allcycles(graph);
% p = plot(graph);
max_cycles_to_plot = 0;
colors = {'red', 'blue', 'green', 'black', 'cyan', 'magenta'};
lineStyles = {"-"}; % {":", "--", "-."};
% return

% ⋘────────── Get the nerighbors of a vertex ──────────⋙
v0_ndx = 1;
v0 = conical_partition.getStateSpaceVertex(v0_ndx);
v0_nbd_ndxs = conical_partition.getVerticesAdjacentToVertex(v0_ndx);
adjacent_regions_ndxs = conical_partition.getAdjacentRegionsIndices(v0_ndx);
v_prev_ndx = v0_nbd_ndxs(1);
v_next_ndx = v0_nbd_ndxs(2);
v_prev = conical_partition.getStateSpaceVertex(v_prev_ndx);
v_next = conical_partition.getStateSpaceVertex(v_next_ndx);
v0v_prev_ndxs = conical_partition.getConjoiningRegionsIndices(v0_ndx, v_prev_ndx);
v0v_next_ndxs = conical_partition.getConjoiningRegionsIndices(v0_ndx, v_next_ndx);

% ⋘────────── Check if v_prev is reachable from v0 ──────────⋙
D_prev = conical_partition.getDerivativeCone(v0v_prev_ndxs);
C_prev = conical_partition.getStateCone(v0v_prev_ndxs);
R_prev = (v0 + 1e3 * D_prev);
% R_prev.removeVertex(v0);
% disp("Intersection (R_" + v0v_prev_ndxs + " ∩ v_" + v_prev_ndx + "): ")
% disp(R_prev.intersectRay(v_prev))

% ⋘────────── Check if v_next is reachable from v0 ──────────⋙
D_next = conical_partition.getDerivativeCone(v0v_next_ndxs);
C_next = conical_partition.getStateCone(v0v_next_ndxs);
R_next = (v0 + 1e3 * D_next);
% disp('is member')
% ismember(v0', R_next.vertices', 'rows')
% disp('is member with tolerance')
% [~, v0ndx] = ismembertol(v0', R_next.vertices', 'ByRows', true);

% disp("Intersection (R_" + v0v_next_ndxs + " ∩ v_" + v_next_ndx + "):")
% disp(R_next.intersectRay(v_next))
% v_left = conical_partition.getStateSpaceVertex(v0_nbd_ndxs(1))
% v_next = conical_partition.getStateSpaceVertex(v0_nbd_ndxs(2))
% D.linearTransform(1e4 * A_cont)
% % (v0 + (1e4 * A_cont) * D) | C
% % (v0 + (1e4 * A_cont) * D) | [[0; 0], 1000*v_left]
% R = (v0 + (1e1 * A_cont) * D);

% D = conical_partition.getDerivativeCone(v0_ndx);

% D
% (v0 + D.linearTransform(1e4 * A_cont))
% angle_array = conical_partition.stateVertexAngles;

% cone_angle = angle_array(2) - angle_array(1);

% ╭────────────────────────────────────╮
% │             Plot Graph             │
% ╰────────────────────────────────────╯
pwintz.plots.namedFigure("Graph");
clf;
hold on;
xlim(1.5*[-1, 1]);
ylim(1.5*[-1, 1]);
axis square;
pwintz.plots.plotUnitCircle();

graph = conic_abstraction.graph;
p = plot(graph, ...
  "XData", graph.Nodes.x, ...
  "YData", graph.Nodes.y, ...
  "NodeColor", "black", ...
  "EdgeAlpha", 1, ...
  'NodeLabel', graph.Nodes.TeXLabel, ...
  'EdgeLabel', graph.Edges.WeightIntervalStr, ...
  "NodeFontSize", 12, ...
  "EdgeFontSize", 10 * graph.Edges.Length / max(graph.Edges.Length), ...
  "Interpreter", "latex"...
);
% ⋘────────── Visualize the weights ──────────⋙
geq_1_color = [1 0 0];
lt_1_color = [0 0 1];
edge_color = (graph.Edges.Weight >= 1) * geq_1_color + (graph.Edges.Weight < 1) * lt_1_color;
weight_range = max(abs(graph.Edges.Weight - 1));
% p.LineWidth = 5 * (1 + (graph.Edges.Weight - 1) / weight_range);
% p.LineWidth = 5*graph.Edges.Weight / weight_range;

% For the plotted weights, we want to illustrate that the importance aspect of an edges weight is its distance above or below 1. We use abs(log(weight)) to normalize. This way, weights equal to 1/2 and 2 are plotted with the same width.
weight = abs(log(graph.Edges.Weight));
weight = 4 * weight / max(weight) + 1;
p.LineWidth = weight;
p.ArrowSize = 8*sqrt(p.LineWidth);
p.EdgeColor = edge_color;

% ⋘────────── Highlight cycles ──────────⋙
for i = 1:min(numel(cycles), max_cycles_to_plot)
  color = colors{mod(i, length(colors)) + 1};
  lineStyle = lineStyles{mod(i, length(lineStyles)) + 1};
  highlight(p,'Edges', edgecycles{i}, 'EdgeColor', color, 'LineStyle', lineStyle, 'NodeColor', color, 'LineWidth',1.5, 'MarkerSize',6);
end


pwintz.plots.namedFigure("2D Vertices");
clf;
hold on;
xlim(1.5*[-1, 1]);
ylim(1.5*[-1, 1]);
axis square;

drawRay([0 0 v0'],     "Color", "black", 'LineWidth', 2, "DisplayName", "ray$(v_{" + v0_ndx + "})$");
drawRay([0 0 v_prev'], "Color", 0.8*[0 0 1], 'LineWidth', 2, "DisplayName", "ray$(v_{" + v_prev_ndx + "})$");
drawRay([0 0 v_next'], "Color", 0.8*[1 0 0], 'LineWidth', 2, "DisplayName", "ray$(v_{" + v_next_ndx + "})$");
fillPolygon(C_prev.vertices', "FaceColor", 0.6*[0 0 1], "DisplayName", "$C_{" + v_prev_ndx + "}$", "LineWidth", 3);
fillPolygon(C_next.vertices', "FaceColor", 0.6*[1 0 0], "DisplayName", "$C_{" + v_next_ndx + "}$", "LineWidth", 3);
fillPolygon(R_prev.vertices', "FaceColor", 0.9*[0 0 1], "DisplayName", "$R_{" + v_prev_ndx + "}$", "LineWidth", 3);
fillPolygon(R_next.vertices', "FaceColor", 0.9*[1 0 0], "DisplayName", "$R_{" + v_next_ndx + "}$", "LineWidth", 3);

R_prev_cap_C_prev = R_prev | 1e4*C_prev;
R_next_cap_C_next = R_next | 1e4*C_next;
fill_prev_args = {"FaceColor", 0.7*[0 0 1], "FaceAlpha", 1, "DisplayName", "$R_{" + v_prev_ndx + "} \cap C_{" + v0v_prev_ndxs + "}$", "LineWidth", 3};
fill_next_args = {"FaceColor", 0.7*[1 0 0], "FaceAlpha", 1, "DisplayName", "$R_{" + v_next_ndx + "} \cap C_{" + v0v_next_ndxs + "}$", "LineWidth", 3};
if ~isempty(R_prev_cap_C_prev)
  fillPolygon(R_prev_cap_C_prev.vertices', fill_prev_args{:});
  % R_prev_cap_ray_v_prev = R_prev_cap_C_prev.removeVertex(v0)
  % vecnorm(R_prev_cap_ray_v_prev.vertices)% norm of each column
else
  % TODO: Add an invisible plot that appears in the legend.
  % plot([nan nan nan], [nan nan nan], fill_prev_args{:})
end
if ~isempty(R_next_cap_C_next)
  fillPolygon(R_next_cap_C_next.vertices', fill_next_args{:});
  %  R_next_cap_ray_v_next = R_next_cap_C_next.removeVertex(v0)
  % vecnorm(R_next_cap_ray_v_next.vertices)% norm of each column
else
  % TODO: Add an invisible plot that appears in the legend.
  % plot([nan nan nan], [nan nan nan], fill_next_args{:})
end

% fillPolygon(R_next.vertices' | , "DisplayName", "$R_{" + v_next_ndx + "}$", "LineWidth", 3);
% fillPolygon(R_next.vertices' | , "DisplayName", "$R_{" + v_next_ndx + "}$", "LineWidth", 3);


legend("interpreter", "latex", "fontSize", 12);

% ╭──────────────────────────────────────────────────╮
% │             Plot Linear Vector Field             │
% ╰──────────────────────────────────────────────────╯
pwintz.plot.plotLinearVectorField(A_cont, DisplayName="$x \mapsto Ax$");

% ⋘────────── Plot circle ──────────⋙
pwintz.plots.plotUnitCircle();

graph = conic_abstraction.graph;
p = plot(graph, "XData", graph.Nodes.x, "YData", graph.Nodes.y);

for i = 1:numel(cycles)
  color = colors{mod(i, length(colors)) + 1};
  highlight(p,'Edges', edgecycles{i}, 'EdgeColor', color, 'NodeColor', color, 'LineWidth',1.5, 'MarkerSize', 6);
end

return

% % ⋘────────── Inscribing Polygon ──────────⋙
% plot(cos(angle_array), sin(angle_array), '-r', DisplayName='Inscribing Polygon');
% 
% % ⋘────────── Circumscribing Polygon ──────────⋙
% r_circumcircle = sec(conical_partition.stateVertexAngles / 2);
% plot(r_circumcircle.*cos(angle_array), r_circumcircle.*sin(angle_array), '-b', DisplayName='Circumscribing Polygon');
% legend("AutoUpdate", false)
% clear theta

r = 1;
for i = 1:conical_partition.nSlices
  % x = cos(theta);
  % y = sin(theta);
  % vDeriv = conical_partition.getDerivativeSpaceVertex(i);
  vState = conical_partition.getStateSpaceVertex(i);
  plot(r*[0, vState(1)], r*[0, vState(2)], "Marker", "none", "LineStyle", "-.", "Color", "black");
end

for i = 1:conical_partition.nSlices
  v1 = conical_partition.getStateSpaceVertex(i);
  v2 = conical_partition.getStateSpaceVertex(i+1);
  angle = conical_partition.getStateSpaceVertexAngle(i+1) - conical_partition.getStateSpaceVertexAngle(i);

  % We scale the outer edge of the polynomial by the "circumscribe_radius", which is the smallest r > 0 such that if we multiply each vertex by r, then the line between the resulting points does only touches the boundary of the unit disk. 
  % Given theta in (-pi/2, pi/2) as the angle between the vertices, the circumscribe_radius is hypotenuse of a right triangle formed as follows:
  % 1. Consider the triangle formed by v1, 0, v2. 
  % 2. Let b be the point where the bisection of v1, 0, v2 intersects the unit circle. 
  % 3. Then, v1, 0, b is a right triangle with 
  %   * Hypotenuse equal to r="circumscribe_radius"
  %   * A_cont side with length 1 that goes from the origin to b.
  %   * The adjacent angle is theta/2. 
  % Thus, cos(theta/2) = (hypotenuse)/(adjacent) = 1 / r. 
  % 4. Solving for r produces  r = sec(theta/2).
  circumscribe_radius = abs(sec(angle / 2));
  P = [v1, circumscribe_radius * v1, circumscribe_radius * v2, v2];
  plotPolygon(P)

  C = conical_partition.getStateCone(i);
  plotPolygon(C, color="blue", alpha=0.1);

  D = conical_partition.getDerivativeCone(i);
  plotPolygon(C, color="blue", alpha=0.1);

  % Compute the reachable set from P flowing according to \dot x \in AC = D.
  % R = addDerivativeCone(P, 1000*D);
  % break
  % pause(0.01)

  % R_in_C = intersection(R, 1000*C);
  % plotPolygon(R_in_C, color="green", alpha=0.1)
  
  % ╭──────────────────────────────────────────────────╮
  % │             Reachability from Vertex             │
  % ╰──────────────────────────────────────────────────╯
  R = addDerivativeConeToStateVertex(v1, 1000*D);
  R_in_C = intersection(R, 1000*C);
  plotPolygon(R_in_C, color=[0.9 0.0 0.3], alpha=0.8)
  R = addDerivativeConeToStateVertex(v2, 1000*D);
  R_in_C = intersection(R, 1000*C);
  plotPolygon(R_in_C, color=[0.9 0.3 0], alpha=0.8)
end

% function result = linearTransform(A_cont, vertices)
%   result = 
% end

function result = addDerivativeConeToStateVertex(v, D)
  % Given a convex polygon P, add the derivative cone D, which consists of three points [v1, 0, v2]. 
  % The addition is in the sense of the Minkowski sum.
  v1 = D(:, 1);
  v2 = D(:, 3);
  result = v + D;
  % dTriangulated = delaunayTriangulation(result')
  % convexHullIndices = convexHull(dTriangulated)
  % result = dTriangulated.Points(convexHullIndices,:)'
  % % Use geom3d trim extra points
  % result = convexHull(result);
end

function plotPolygon(vertices, varargin)
  % Parse inputs from varargin argument.
  p = inputParser;
  addParameter(p, "Color", "r"); % Name-value parameter.
  addParameter(p, "Alpha", 0.5); % Name-value parameter.
  parse(p, varargin{:});
  args = p.Results;
  patch(vertices(1,:), vertices(2,:), args.Color, "FaceAlpha", args.Alpha)
end

