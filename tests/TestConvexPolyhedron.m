classdef TestConvexPolyhedron < matlab.unittest.TestCase
  % ! You can run these tests with 
  % ` runtests TestConvexPolyhedron

  % ╭───────────────────────────────────╮
  % │ ╭───────────────────────────────╮ │
  % │ │             Tests             │ │
  % │ ╰───────────────────────────────╯ │
  % ╰───────────────────────────────────╯
  methods (Test)
    % ╭───────────────────────────────────────────────────╮
    % │             Equality "==" overloading             │
    % ╰───────────────────────────────────────────────────╯
    function test_eq_basicEquality(testCase)
      % ⋘────────── Setup ───────────⋙
      cp1_equal    = Polytope.fromConvexHull([[1; 1], [0; 0], [-1; 2]]);
      cp2_equal    = Polytope.fromConvexHull([[1; 1], [0; 0], [-1; 2]]);
      cp3_notEqual = Polytope.fromConvexHull([[1; 1], [0; 3], [-1; 2]]);
      
      % ⋘────────── Execute and  Verify ──────────⋙
      testCase.assertEqual(cp1_equal, cp1_equal, "same object");
      testCase.assertEqual(cp1_equal, cp2_equal);
      testCase.assertNotEqual(cp1_equal, cp3_notEqual);
      % Test operators correctly produce "true".
      testCase.assertTrue( cp1_equal == cp1_equal,    "Same object");
      testCase.assertTrue( cp1_equal == cp2_equal,    "Different object, but equal");
      testCase.assertTrue( cp1_equal ~= cp3_notEqual, "Different objects, not equal");
      % Test operators correctly produce "false".
      testCase.assertFalse(cp1_equal == cp3_notEqual, "Not equal with ""=="" operator");
      testCase.assertFalse(cp1_equal ~= cp2_equal,    "Equal with ""~="" operator");
    end % End of function.

    function test_eq_permutedVertexOrder(testCase)
      % ⋘────────── Setup ───────────⋙
      cp1_equal    = Polytope.fromConvexHull([[1; 1], [0; 0], [-1; 2]]);
      cp2_equal    = Polytope.fromConvexHull([[0; 0], [1; 1], [-1; 2]]);
      cp3_notEqual = Polytope.fromConvexHull([[1; 1], [0; 3], [-1; 2]]);
      
      % ⋘────────── Execute and  Verify ──────────⋙
      % Test operators correctly produce "true".
      testCase.assertTrue( cp1_equal == cp1_equal,    """=="" operator is reflexive");
      testCase.assertTrue( cp1_equal == cp2_equal,    pwintz.strings.format("Equivalent objects %s and %s are equal with ""=="" operator", cp1_equal, cp2_equal));
      testCase.assertTrue( cp1_equal ~= cp3_notEqual, "Not equal with ""~="" operator");
      % Test operators correctly produce "false".
      testCase.assertFalse(cp1_equal == cp3_notEqual, "Not equal with ""=="" operator");
      testCase.assertFalse(cp1_equal ~= cp1_equal,    "Equal with ""=="" operator");
    end % End of function.

    % ╭───────────────────────────────────────────╮
    % │             "contains" method             │
    % ╰───────────────────────────────────────────╯
    function test_contains_singlePoint(testCase)
      % ⋘────────── Setup ───────────⋙
      p = Polytope.fromConvexHull([[10; 0], [ 0; 0], [0; 10], [10; 10]]);

      % ⋘────────── Execute and Verify ──────────⋙
      testCase.verifyTrue(p.containsPoints([1; 1]),   "Interior point");
      testCase.verifyTrue(p.containsPoints([10; 0]),  "Vertex");
      testCase.verifyFalse(p.containsPoints([-1; 0]), "Exterior point");
      
      % ! Testing the inclusion of boundary points may be unreilable due to numerical issues, but it should work in this case. 
      testCase.verifyTrue(p.containsPoints([1; 0]), "Boundary point (not vertex)");
    end % End of function.


    function test_contains_multiplePoints(testCase)
      % ⋘────────── Setup ───────────⋙
      p = Polytope.fromConvexHull([[10; 0], [ 0; 0], [0; 10], [10; 10]]);

      % ⋘────────── Execute and Verify ──────────⋙
      testCase.verifyEqual(p.containsPoints([[1; 1], [10; 0], [-1; 0]]), logical([1 1 0]));
    end % End of function.


    % ╭──────────────────────────────────────────────────╮
    % │  ╭────────────────────────────────────────────╮  │
    % │  │             Inersection Method             │  │
    % │  ╰────────────────────────────────────────────╯  │
    % ╰──────────────────────────────────────────────────╯
    function test_intersection_polyhedron_cap_polyhedron(testCase)
      % ⋘────────── Setup ───────────⋙
      h1 = HalfspaceRepresentation([-1,  0], [0]);
      h2 = HalfspaceRepresentation([ 0, -1], [0]);

      p1 = ConvexPolyhedron(h1);
      p2 = ConvexPolyhedron(h2);
      
      % ⋘────────── Execute ─────────⋙
      
      p_cap_p = intersection(p1, p2);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertInstanceOf(p_cap_p, "ConvexPolyhedron");
      testCase.assertTrue(p_cap_p.containsPoints([1; 1]));
      testCase.assertFalse(p_cap_p.containsPoints([-1; 1]));
      testCase.assertFalse(p_cap_p.containsPoints([1; -1]));
    end % End of function.

    function test_intersection_Polyhedron_cap_PolyhedralCone_(testCase)
      % ⋘────────── Setup ───────────⋙
      A_ineq_1 =[-0.4, -0.92; 
                  0.29, 0.96];
      b_ineq_1 =[0.4; -0.29];
      h1 = HalfspaceRepresentation(A_ineq_1, b_ineq_1);
      cp = ConvexPolyhedron(h1);

      rays = [
        -1.0000   -0.9980; 
        0.0000     0.0628;
      ];
      cone = ConvexPolyhedralCone.fromRays(rays);

      % h2 = cone.halfspace_representation;

      h1
      figure(1);
      clf();
      xlim(3*[-1, 1]);
      ylim(3*[-1, 1]);
      axis square;
      hold on;
      h1.plot();
      cone.plot();

      % ⋘────────── Execute ─────────⋙
      
      p_cap_p = intersection(cp, cone);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertInstanceOf(p_cap_p, "ConvexPolyhedron");
      testCase.assertTrue(p_cap_p.containsPoints([1; 1]));
      testCase.assertFalse(p_cap_p.containsPoints([-1; 1]));
      testCase.assertFalse(p_cap_p.containsPoints([1; -1]));
    end % End of function.

    % ╭─────────────────────────────────────────────╮
    % │             atXDoesVPointInward             │
    % ╰─────────────────────────────────────────────╯
    function test_atXDoesVPointInward_pointOnBoundary_2D(testCase)
      % ⋘────────── Setup ───────────⋙
      % A = -eye(2); % First quandrant.
      % b = [0; 0];
      rays = [[1; 0], [0; 1]];
      poly = ConvexPolyhedralCone.fromRays(rays);
      
      % ⋘────────── Execute and Verify ──────────⋙
      testCase.assertTrue(poly.atXDoesVPointInward([0; 0], [1; 1]));
      testCase.assertTrue(poly.atXDoesVPointInward([1; 0], [-1; 1]));
      testCase.assertTrue(poly.atXDoesVPointInward([0; 1], [1; -2]));
      
      testCase.assertFalse(poly.atXDoesVPointInward([0; 0], [-1; -1]));
      testCase.assertFalse(poly.atXDoesVPointInward([0; 0], [1; -1]));
      testCase.assertFalse(poly.atXDoesVPointInward([1; 0], [0; -1]));
      testCase.assertFalse(poly.atXDoesVPointInward([0; 1], [-1; 0]));
    end % End of function.

  end % End of test methods

end
