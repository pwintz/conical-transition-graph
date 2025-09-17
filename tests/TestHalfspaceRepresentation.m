classdef TestHalfspaceRepresentation < matlab.unittest.TestCase
% ` runtests TestHalfspaceRepresentation
  % ╭───────────────────────────────────╮
  % │ ╭───────────────────────────────╮ │
  % │ │             Tests             │ │
  % │ ╰───────────────────────────────╯ │
  % ╰───────────────────────────────────╯
  methods (Test)
    function test_constructor_empty(testCase)
      % ⋘────────── Execute ─────────⋙
      h = HalfspaceRepresentation();
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(h.n_inequality_constraints, int32(0));
      testCase.assertEqual(h.n_equality_constraints,   int32(0));
      testCase.assertFalse(h.isBounded());
      
    end % End of function.
    
    function test_constructor_noEqualityConstraints(testCase)
      % ⋘──────── Setup ────────⋙
      rng(1);
      n_ineq = 2;
      n_eq   = 0;
      dim    = 2;
      A_ineq = rand(n_ineq, dim);
      b_ineq = ones(n_ineq, 1);
    
      % ⋘────────── Execute ─────────⋙
      h = HalfspaceRepresentation(A_ineq, b_ineq);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(h.n_inequality_constraints, int32(n_ineq));
      testCase.assertEqual(h.n_equality_constraints,   int32(n_eq));
      testCase.assertEqual(h.ambient_dimension,        int32(dim));
      testCase.assertFalse(h.isBounded());
      
    end % End of function.
    
    function test_constructor_withEqualityConstraints(testCase)
      % ⋘──────── Setup ────────⋙
      n_ineq = 2;
      n_eq   = 5;
      dim    = 6;
      rng(1); % Set seed for random number generator.
      A_ineq = rand(n_ineq, dim);
      b_ineq = rand(n_ineq, 1);
      A_eq = rand(n_eq, dim);
      b_eq = rand(n_eq, 1);
    
      % ⋘────────── Execute ─────────⋙
      h = HalfspaceRepresentation(A_ineq, b_ineq, A_eq, b_eq);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(h.n_inequality_constraints, int32(n_ineq));
      testCase.assertEqual(h.n_equality_constraints,   int32(n_eq));
      testCase.assertEqual(h.ambient_dimension,        int32(dim));
      testCase.assertFalse(h.isBounded());
      
    end % End of function.
    
    
    % ╭────────────────────────────────────────────────────────╮
    % │             Constructor Merges Constraints             │
    % ╰────────────────────────────────────────────────────────╯
    function test_constructor_normalizes_A_and_b(testCase)
      % ⋘────────── Setup ───────────⋙
      A_ineq = [1, 2];
      b_ineq = 3;
      A_eq = [4, 5];
      b_eq =    -2;
      
      % ⋘────────── Execute ─────────⋙
      h = HalfspaceRepresentation(A_ineq, b_ineq, A_eq, b_eq);
      
      % ⋘────────── Verify ──────────⋙
      A_ineq_expected = [1, 2] / norm([1, 2]);
      b_ineq_expected = 3 / norm([1, 2]);
      testCase.assertEqual(h.n_inequality_constraints, int32(1));
      testCase.assertEqual(h.A_ineq, A_ineq_expected, "A_ineq");
      testCase.assertEqual(h.b_ineq, b_ineq_expected, "b_ineq");
      
    end % End of function.

    function test_constructor_removesIdenticalInequalityConstraints(testCase)
      % ⋘────────── Setup ───────────⋙
      A_ineq = [1, 2; 1, 2];
      b_ineq = [3; 3];
      
      % ⋘────────── Execute ─────────⋙
      h = HalfspaceRepresentation(A_ineq, b_ineq);
      
      % ⋘────────── Verify ──────────⋙
      A_ineq_expected = [1, 2] / norm([1, 2]);
      b_ineq_expected = 3 / norm([1, 2]);
      testCase.assertEqual(h.n_inequality_constraints, int32(1));
      testCase.assertEqual(h.A_ineq, A_ineq_expected, "A_ineq");
      testCase.assertEqual(h.b_ineq, b_ineq_expected, "b_ineq");
      
    end % End of function.

    % ╭──────────────────────────────────╮
    % │      Constructor Exceptions      │
    % ╰──────────────────────────────────╯
    function test_constructor_zeroRowCausesError(testCase)

      A_ineq = [1, 2; 
                0, 0]; % <- Zero row.
      b_ineq = [3; 3];
      
      % ⋘────────── Execute ─────────⋙
      testCase.assertError(@() HalfspaceRepresentation(A_ineq, b_ineq), "HalfspaceRepresentation:mustHaveOnlyNonzeroRows");
      
    end % End of function.

    % ╭─────────────────────────────────────╮
    % │             isBounded()             │
    % ╰─────────────────────────────────────╯
    function test_isBounded_1D_bounded(testCase)
      % ⋘──────── Setup ────────⋙
      n_ineq = 2;
      dim    = 1;
      A_ineq = [ 1; 
                -1];
      b_ineq = [ 1;
                 0];
    
      % Check that the constraints define the interval [0, 1]
      assert(all(A_ineq * [ 0] <= b_ineq)); %#ok<NBRAK1>
      assert(all(A_ineq * [ 1] <= b_ineq)); %#ok<NBRAK1>
      assert(any(A_ineq * [-1] >= b_ineq)); %#ok<NBRAK1>
      assert(any(A_ineq * [ 2] >= b_ineq)); %#ok<NBRAK1>

      % ⋘────────── Execute ─────────⋙
      h = HalfspaceRepresentation(A_ineq, b_ineq);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(h.n_inequality_constraints, int32(n_ineq));
      testCase.assertEqual(h.ambient_dimension,        int32(dim));
      testCase.assertTrue(h.isBounded());
      
    end % End of function.

    function test_isBounded_4D_unbounded(testCase)
      % ⋘──────── Setup ────────⋙
      n_ineq = 5;
      dim    = 4;
      A_ineq = ones(n_ineq, dim);
      b_ineq = ones(n_ineq, 1);
      h = HalfspaceRepresentation(A_ineq, b_ineq);
    
      % ⋘────────── Execute and Verify ──────────⋙
      testCase.assertFalse(h.isBounded());
      
    end % End of function.
    
    function test_intersection(testCase)
      % ⋘──────── Setup ────────⋙
      rng(1);
      dim     = 7;
      n_ineq1  = 3;
      n_ineq2  = 2;
      A1_ineq = randn(n_ineq1, dim);
      b1_ineq = randn(n_ineq1, 1);
      A2_ineq = randn(n_ineq2, dim);
      b2_ineq = randn(n_ineq2, 1);
      h1 = HalfspaceRepresentation(A1_ineq, b1_ineq);
      h2 = HalfspaceRepresentation(A2_ineq, b2_ineq);
    
      % ⋘────────── Execute ─────────⋙
      h_inter = h1.intersection(h2);
      
      % ⋘────────── Verify ──────────⋙
      n_tests = 100;
      for x = randn(dim, n_tests)
        % Randomly sample points, then test whether they are in the intersection, explicitly.
        is_in_h1      = h1.containsPoint(x);
        is_in_h2      = h2.containsPoint(x);
        is_in_h_inter = h_inter.containsPoint(x);
        testCase.assertEqual(is_in_h_inter, is_in_h1 && is_in_h2);
      end

    end % End of function.

    function test_intersection_cone2D(testCase)
      % ⋘──────── Setup ────────⋙
      rng(1);
      dim     = 2;
      A1_ineq = [0, -1];
      b1_ineq = 0;
      A2_ineq = [-1, 0];
      b2_ineq = 0;
      h1 = HalfspaceRepresentation(A1_ineq, b1_ineq);
      h2 = HalfspaceRepresentation(A2_ineq, b2_ineq);
    
      % ⋘────────── Execute ─────────⋙
      h_inter = h1.intersection(h2);
      
      % ⋘────────── Verify ──────────⋙
      n_tests = 100;
      for x = randn(dim, n_tests)
        % Randomly sample points, then test whether they are in the intersection, explicitly.
        is_in_h1      = h1.containsPoint(x);
        is_in_h2      = h2.containsPoint(x);
        is_in_h_inter = h_inter.containsPoint(x);
        testCase.assertEqual(is_in_h_inter, is_in_h1 && is_in_h2);
      end

      % Check that we cna
      h_inter.getConeRays();
    end % End of function.

    function test_intersection_unbounded3D(testCase)
      % ⋘──────── Setup ────────⋙
      rng(1);
      dim     = 3;
      A1_ineq = [0, -1, 0; 1, -1, 0];
      b1_ineq = [0, 0];
      A2_ineq = [-1, 0, 0];
      b2_ineq = 0;
      h1 = HalfspaceRepresentation(A1_ineq, b1_ineq);
      h2 = HalfspaceRepresentation(A2_ineq, b2_ineq);
    
      % ⋘────────── Execute ─────────⋙
      h_inter = h1.intersection(h2);
      
      % ⋘────────── Verify ──────────⋙
      n_tests = 100;
      for x = randn(dim, n_tests)
        % Randomly sample points, then test whether they are in the intersection, explicitly.
        is_in_h1      = h1.containsPoint(x);
        is_in_h2      = h2.containsPoint(x);
        is_in_h_inter = h_inter.containsPoint(x);
        testCase.assertEqual(is_in_h_inter, is_in_h1 && is_in_h2);
      end
    end % End of function.

    function test_intersection_unboundedReducedDimension(testCase)
      % ⋘──────── Setup ────────⋙
      rng(1);
      dim     = 3;
      A1_ineq = [0, -1, 0];
      b1_ineq = 0;
      A2_ineq = [0, 1, 0];
      b2_ineq = 0;
      h1 = HalfspaceRepresentation(A1_ineq, b1_ineq);
      h2 = HalfspaceRepresentation(A2_ineq, b2_ineq);
    
      % ⋘────────── Execute ─────────⋙
      h_inter = h1.intersection(h2);
      
      % ⋘────────── Verify ──────────⋙
      n_tests = 100;
      for x = randn(dim, n_tests)
        % Randomly sample points, then test whether they are in the intersection, explicitly.
        is_in_h1      = h1.containsPoint(x);
        is_in_h2      = h2.containsPoint(x);
        is_in_h_inter = h_inter.containsPoint(x);
        testCase.assertEqual(is_in_h_inter, is_in_h1 && is_in_h2);
      end
    end % End of function.

    function test_getConeRays_ifNoConstraints(testCase)
      % ⋘────────── Setup ───────────⋙
      unconstrained_hs = HalfspaceRepresentation();
      
      % ⋘────────── Execute and Verify─────────⋙
      testCase.assertEmpty(unconstrained_hs.getConeRays());
      
    end % End of function.

    function test_getConeRays(testCase)
      % ⋘────────── Setup ───────────⋙
      A = -[-1, 1; 
             1, 0];
      b = [0, 0];
      hs = HalfspaceRepresentation(A, b);
      
      % ⋘────────── Execute ─────────⋙
      rays = hs.getConeRays();
      testCase.assertTrue(pwintz.arrays.isColumnIn([0; 1], rays, tolerance=1e-12));
      testCase.assertTrue(pwintz.arrays.isColumnIn(1/sqrt(2) * [1; 1], rays, tolerance=1e-12));

      % ⋘────────── Verify ─────────⋙
      
    end % End of function.


    function test_getConeRays_withNearlyZero_b_values(testCase)
      % ⋘────────── Setup ───────────⋙
      A = -[-1, 1; 
             1, 0];
      b = [0, 1e-17];
      hs = HalfspaceRepresentation(A, b);
      
      % ⋘────────── Execute ─────────⋙
      rays = hs.getConeRays();
      testCase.assertTrue(pwintz.arrays.isColumnIn([0; 1], rays, tolerance=1e-12));
      testCase.assertTrue(pwintz.arrays.isColumnIn(1/sqrt(2) * [1; 1], rays, tolerance=1e-12));
    
      % ⋘────────── Verify ─────────⋙
      
    end % End of function.

    % ╭───────────────────────────────────────────────╮
    % │             getBoxBoundedVertices             │
    % ╰───────────────────────────────────────────────╯
    function test_getConeRays_degenerateCase(testCase)
      % ⋘────────── Setup ───────────⋙
      A_ineq=[
              -0.17,      -0.056,         0.98;
              0.17,       -0.056,       -0.98;
              -1.4e-16,   -1,       0;
              1.2e-16,     1,       0;
              0.15,        0.47,       0.87;
              0.46, -0.39, -0.8;
            ];
      b_ineq = [-2.8e-17;
                2.4e-17;
                2.2e-18;
                -5.6e-17;
                1.8e-17;
                5.6e-17];
      % b_ineq = 0* b_ineq
      % b_ineq = zeros(size(A_ineq, 1), 1);
      h = HalfspaceRepresentation(A_ineq, b_ineq);
      
      % ⋘────────── Execute ─────────⋙
      h.getConeRays();
      
    end % End of function.

    % ╭──────────────────────────────────────────────╮
    % │  ╭────────────────────────────────────────╮  │
    % │  │             ContainsPoints             │  │
    % │  ╰────────────────────────────────────────╯  │
    % ╰──────────────────────────────────────────────╯
    function test_containsPoints_multiple(testCase)
      % ⋘────────── Setup ───────────⋙
      dimension=2;
      h = HalfspaceRepresentation.getSquare(dimension);
      
      % ⋘────────── Execute ─────────⋙
      points = [[1; 0], [0; 0], [20; 0]];
      is_in_box = h.containsPoints(points);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(is_in_box, logical([1, 1, 0]));
      
      % Check result agains t "containsPoint" method.
      expected = [h.containsPoint(points(:, 1)), h.containsPoint(points(:, 2)), h.containsPoint(points(:, 3))];
      testCase.assertEqual(is_in_box, expected);
      
    end % End of function.


    % ╭───────────────────────────────────╮
    % │             FromConvexHull             │
    % ╰───────────────────────────────────╯
    function test_fromConvexHull_single_point(testCase)
      % ⋘────────── Setup ───────────⋙
      point = [1; 2];

      % ⋘────────── Execute ─────────⋙
      h = HalfspaceRepresentation.fromConvexHull(point);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertTrue(h.containsPoint(point));
    end % End of function.

    function test_fromConvexHull_caseThatWasProducingErrors(testCase)
      % ⋘────────── Setup ───────────⋙
      points = [ 
        1.0000,    0.9980,    1.0005,    0.9985;
             0,    0.0628,         0,    0.0628;
      ];
      pwintz.plots.plotPoints(points, "plotArgs", {"LineStyle", "-"});
      
      % ⋘────────── Execute ─────────⋙
      h = HalfspaceRepresentation.fromConvexHull(points);
      bounded_vertices = h.getBoxBoundedVertices();
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertTrue(h.isBounded());
    end % End of function.

    function test_fromConvexHull_anotherCausingProblems(testCase)
      % ⋘────────── Setup ───────────⋙
      points = [
        1.0000,  0.9980,  1.0005,  0.9985,  0.3249,  0.3229,  0.3254,  0.3234,  0.3371,  0.3352,  0.3376,  0.3357;
        0, 0.0628, 0, 0.0628   -0.7377   -0.6749   -0.7377   -0.6749   -0.7487   -0.6860   -0.7487   -0.6859;
      ];
      
      % ⋘────────── Execute ─────────⋙
      h = HalfspaceRepresentation.fromConvexHull(points);
      bounded_vertices = h.getBoxBoundedVertices();
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertTrue(h.isBounded());

    end % End of function.


    % ╭─────────────────────────────────────────╮
    % │             fromConicalHull             │
    % ╰─────────────────────────────────────────╯
    function test_fromConicalHull_nearlyParallelRays(testCase)
      % ⋘────────── Setup ───────────⋙
      rays = [
        -0.98228725, -0.98006658; 
         0.18738131,  0.19866933; 
      ];
      % rays(:, 1) - rays(:, 2)
      rays = [
        [1; 1],  [1.01; 1];
      ];
      
      % ⋘────────── Execute ─────────⋙
      
      hs = HalfspaceRepresentation.fromConicalHull(rays);
      
      % ⋘────────── Verify ──────────⋙
      % testCase.verifyFail("Test case needs to be implemented.");
    end % End of function.

    function test_fromConicalHull_2_parallel_rays(testCase)
      % ⋘────────── Setup ───────────⋙
      ray = [1; 1];
      rays = [ray, 2*ray];
      
      % ⋘────────── Execute ─────────⋙
      hs = HalfspaceRepresentation.fromConicalHull(rays);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertTrue(hs.containsPoint(ray));
      testCase.assertTrue(hs.containsPoints(ray));

      % The negation of the ray should not be in the set.
      testCase.assertFalse(hs.containsPoint(-ray));
      
      % Only one ray is kept
      testCase.assertEqual(hs.rays, ray ./ norm(ray)); 
    end % End of function.

    % ╭────────────────────────────────────────────────────────────╮
    % │             activeInequalityConstraintIndices             │
    % ╰────────────────────────────────────────────────────────────╯
    function test_activeInequalityConstraintIndices(testCase)
      % ⋘────────── Setup ───────────⋙
      % A = -eye(2); % First quandrant.
      % b = [0; 0];
      rays = [[1; 0], [0; 1]];
      h = HalfspaceRepresentation.fromConicalHull(rays);

      horizontal_constraint_ndx = pwintz.arrays.findRowIn([-1, 0], h.A_ineq);
      vertical_constraint_ndx = pwintz.arrays.findRowIn([0, -1], h.A_ineq);
      
      % ⋘────────── Execute and Verify ──────────⋙
      testCase.assertEqual(h.activeInequalityConstraintIndices([0; 0]), [1; 2]);
      testCase.assertEqual(h.activeInequalityConstraintIndices([1; 0]), vertical_constraint_ndx);
      testCase.assertEqual(h.activeInequalityConstraintIndices([0; 1]), horizontal_constraint_ndx);
      testCase.assertEqual(h.activeInequalityConstraintIndices([1; 1]), double.empty(0, 1));
    end % End of function.

    % ╭────────────────────────────────────────────────╮
    % │             String representations             │
    % ╰────────────────────────────────────────────────╯
    function test_char(~)
      % ⋘────────── Setup ───────────⋙
      rays = [[1; 0], [0; 1]];
      
      % ⋘────────── Execute ─────────⋙
      hs = HalfspaceRepresentation.fromConicalHull(rays);
      
      % ⋘────────── Execute ─────────⋙
      char(hs);
      sprintf("%s", hs);
    end % End of function.

  end % End of test methods
end % End of test methods block

