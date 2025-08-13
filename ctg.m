clear

% roll_array = linspace(0, 2*pi, 10)
% pitch_array = linspace(-pi, pi, 10)
% yaw_array = linspace(0, pi, 10);

% ╭────────────────────────────────────────╮
% │ ╭────────────────────────────────────╮ │
% │ │             Given Data             │ │
% │ ╰────────────────────────────────────╯ │
% ╰────────────────────────────────────────╯
% Create a random 2x2 matrix to define \dot x = Ax.
A = rand(2) - 0.5;
% ⋘────────── Some interesting choices of A ──────────⋙
% # This A gives a very large slices in the state space. One of its eigenvalues is -0.007, so the problem likely arises from the nearly zero eigenvalue.
% A =
%     0.4027   -0.0091
%     0.4448   -0.0107

% ╭────────────────────────────────────────────╮
% │ ╭────────────────────────────────────────╮ │
% │ │             CTG Parameters             │ │
% │ ╰────────────────────────────────────────╯ │
% ╰────────────────────────────────────────────╯

conical_partition = ConicalPartition(A, "nSlices", 150);
% angle_array = conical_partition.stateVertexAngles;

% cone_angle = angle_array(2) - angle_array(1);


% pwintz.utils.namedFigure("2D Vertices");
figure(1)
clf
hold on
xlim(1.5*[-1, 1])
ylim(1.5*[-1, 1])
axis square


% ╭──────────────────────────────────────────────────╮
% │             Plot Linear Vector Field             │
% ╰──────────────────────────────────────────────────╯

[x1, x2] = meshgrid(-2:0.1:2,-2:0.1:2);
p = A*[reshape(x1, 1, []); 
       reshape(x2, 1, [])];
px1 = reshape(p(1,:), size(x1));
px2 = reshape(p(2,:), size(x1));
quiver(x1,x2,p1,p2, DisplayName='')

% ⋘────────── Plot circle ──────────⋙
theta = linspace(0, 2*pi, 100);
plot(cos(theta), sin(theta), '-k', DisplayName='Unit Circle')

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
  %   * A side with length 1 that goes from the origin to b.
  %   * The adjacent angle is theta/2. 
  % Thus, cos(theta/2) = (hypotenuse)/(adjacent) = 1 / r. 
  % 4. Solving for r produces  r = sec(theta/2).
  circumscribe_radius = abs(sec(angle / 2));
  P = [v1, circumscribe_radius * v1, circumscribe_radius * v2, v2];
  plotPolygon(P)

  C = conical_partition.getCone(i);
  plotPolygon(C, color="blue", alpha=0.1);

  D = conical_partition.getDerivativeCone(i);
  plotPolygon(C, color="blue", alpha=0.1);

  % Compute the reachable set from P flowing according to \dot x \in AC = D.
  R = addDerivativeCone(P, 1000*D);
  % break
  % pause(0.01)

  R_in_C = intersection(R, 1000*C);
  plotPolygon(R_in_C, color="green", alpha=0.1)
  
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

% function result = linearTransform(A, vertices)
%   result = 
% end

function result = addDerivativeCone(P, D)
  % Given a convex polygon P, add the derivative cone D, which consists of three points [v1, 0, v2]. 
  % The addition is in the sense of the Minkowski sum.
  v1 = D(:, 1);
  v2 = D(:, 3);
  result = [P, P + v1, P + v2];
  dTriangulated = delaunayTriangulation(result');
  convexHullIndices = convexHull(dTriangulated);
  result = dTriangulated.Points(convexHullIndices,:)';
  % Use geom3d trim extra points
  % result = convexHull(result);
end

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

function intersection_vertices = intersection(vertices1, vertices2)
  vertices1
  poly1 = polyshape(vertices1(1,:), vertices1(2,:))
  poly2 = polyshape(vertices2(1,:), vertices2(2,:))
  poly_intersection = intersect(poly1, poly2);
  intersection_vertices = poly_intersection.Vertices';
end