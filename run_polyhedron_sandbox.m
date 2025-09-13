

K = ConvexPolyhedralCone.fromRays([[1; 0; 0], [0; 1; 0], [0; 0; 1]])

vertices = [[1; 2; 2], [0; 1; -1], [3; -1; 2]]
P = Polytope.fromConvexHull(vertices)
assert(P.isBounded())


A = -eye(3);

A*K
% K + K
K + P
P + K

R = intersection(P + A*K, K)


Polytope.getBox(7)
return

pwintz.plots.namedFigure("Reachable set");
plot(K)
clf
plot(R)