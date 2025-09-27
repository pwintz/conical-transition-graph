timeit(@test)

function  test()
  A = rand(2);
  % P1 = ConvexPolyhedron.fromConvexHull([0, 10, 10, 0, 5; 0, 0, 10, 10, 5]);
  % P2 = ConvexPolyhedron.fromConvexHull([0, 2, 10, 0, -1; 0, 0, 10, 2, 15]);
  P = ConvexPolyhedron.fromConvexHull(randn(2, 2000));
  % out = P1 + P2;
  % (A * A) * P;
  A * (A  * (A  * (A  * (A  * (A  * (A  * (A * P)))))));
end