classdef TestConicalPartition < matlab.unittest.TestCase
  % ! You can run these tests using ConicalPartition.test();

%     methods (TestMethodSetup)
%         function setup(testCase)
%             % Setup code
%             rng(1);
%         end
%     end
% 
%     methods (TestMethodTeardown)
%         function teardown(testCase)
%             % Teardown code
%         end
%     end

  % ╭───────────────────────────────────╮
  % │ ╭───────────────────────────────╮ │
  % │ │             Tests             │ │
  % │ ╰───────────────────────────────╯ │
  % ╰───────────────────────────────────╯
  methods (Test)

    % ╭─────────────────────────────────────╮
    % │             Constructor             │
    % ╰─────────────────────────────────────╯

    % ╭─────────────────────────────────────────────────────────────╮
    % │             fromNumberOfUniformDerivativeSlices             │
    % ╰─────────────────────────────────────────────────────────────╯
    function test_fromNumberOfUniformDerivativeSlices(testCase)
        % ⋘────────── Setup ───────────⋙
        A = [0, -1; 1, 0]; % <- rotation by pi/2
        n_slices = 4;
        conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);

        % ⋘────────── Execute and Verify ──────────⋙
        expected_derivative_vertex_angles = [0, pi/2, pi, 3*pi/2];
        expected_state_vertex_angles      = mod([0, pi/2, pi, 3*pi/2] - pi/2, 2*pi);
        testCase.assertEqual(conical_partition.n_vertices, n_slices);

        derivative_angle_error = pwintz.math.angleDistance(conical_partition.derivative_vertex_angles, expected_derivative_vertex_angles);
        state_angle_error = pwintz.math.angleDistance(conical_partition.state_vertex_angles, expected_state_vertex_angles);
        testCase.assertLessThan(derivative_angle_error, 1e-6);
        testCase.assertLessThan(state_angle_error, 1e-6);
        
        testCase.assertEqual(conical_partition.derivative_vertices,     [[1; 0], [0; 1], [-1; 0], [0; -1]], "AbsTol", 1e-6);
        testCase.assertEqual(conical_partition.state_vertices, [[0; -1], [1; 0], [0; 1], [-1; 0]],          "AbsTol", 1e-6);
    end % End test 

    % ╭────────────────────────────────────╮
    % │             fromAngles             │
    % ╰────────────────────────────────────╯
    function test_fromAngles_withoutExtraStateVertexAngles(testCase)
        % ⋘────────── Setup ───────────⋙
        A = magic(2);
        derivative_vertex_angles = [0, 2, 4];
        
        % ⋘────────── Execute ─────────⋙
        conical_partition = ConicalPartition.fromAngles(A, derivative_vertex_angles, []);
        
        
        % ⋘────────── Verify ──────────⋙
        expected_derivative_vertices =     [[cos(0); sin(0)], [cos(2); sin(2)], [cos(4); sin(4)]];
        expected_state_vertices      = A \ [[cos(0); sin(0)], [cos(2); sin(2)], [cos(4); sin(4)]];
        testCase.assertEqual(conical_partition.n_vertices, 3);
        testCase.assertEqual(conical_partition.derivative_vertices, expected_derivative_vertices);
        testCase.assertEqual(conical_partition.state_vertices, expected_state_vertices);
    end % End of function.

    
    function test_fromAngles_withExtraStateVertexAngles_withATheIdentity(testCase)
        % ⋘────────── Setup ───────────⋙
        A = eye(2);
        derivative_vertex_angles  = [0, 2*pi/3, 4*pi/3];
        extra_state_vertex_angles = [pi, pi/4]; % <- Doesn't need to be sorted.

        % ⋘────────── Execute ─────────⋙
        conical_partition = ConicalPartition.fromAngles(A, derivative_vertex_angles, extra_state_vertex_angles);
        
        % ⋘────────── Verify ──────────⋙
        expected_state_angles = [0, pi/4, 2*pi/3, pi, 4*pi/3];
        expected_state_vertices = pwintz.math.angle2UnitVector(expected_state_angles);
        % Because we are using the identity matrix for A, the state and derivative vertices are the same.
        expected_derivative_angles = expected_state_angles;
        expected_derivative_vertices = expected_state_vertices;

        testCase.assertEqual(conical_partition.n_vertices, 5);
        testCase.assertEqual(conical_partition.state_vertex_angles, expected_state_angles, "AbsTol", 1e-6);
        testCase.assertEqual(conical_partition.derivative_vertex_angles, expected_derivative_angles, "AbsTol", 1e-6);
        testCase.assertEqual(conical_partition.derivative_vertices, expected_derivative_vertices, "AbsTol", 1e-6);
        testCase.assertEqual(conical_partition.state_vertices, expected_state_vertices, "AbsTol", 1e-6);
    end % End of function.

    function test_fromAngles_withExtraStateVertexAngles(testCase)
        % ⋘────────── Setup ───────────⋙
        A = magic(2);
        derivative_vertex_angles  = [0, 5, 2, 4, -2];   % <- Not sorted nor in [0, 2*pi).
        extra_state_vertex_angles = [pi, pi/4, -1, 12]; % <- Not sorted nor in [0, 2*pi).
        extra_state_vertices      = pwintz.math.angle2UnitVector(extra_state_vertex_angles);
        extra_derivative_vertices = A * extra_state_vertices;
        extra_derivative_angles   = mod(pwintz.math.atan2(extra_derivative_vertices), 2*pi);

        % ⋘────────── Execute ─────────⋙
        % expected_derivative_vertices =     [[cos(0); sin(0)], [cos(2); sin(2)], [cos(4); sin(4)]];
        % expected_state_vertices      = A \ [[cos(0); sin(0)], [cos(2); sin(2)], [cos(4); sin(4)]];
        
        conical_partition = ConicalPartition.fromAngles(A, derivative_vertex_angles, extra_state_vertex_angles);
        
        % ⋘────────── Verify ──────────⋙
        expected_derivative_angles = sort(mod([derivative_vertex_angles, extra_derivative_angles], 2*pi));
        expected_derivative_vertices = pwintz.math.angle2UnitVector(expected_derivative_angles);
        expected_state_vertices      = A \ expected_derivative_vertices;
        expected_state_angles   = mod(pwintz.math.atan2(expected_state_vertices), 2*pi);

        testCase.assertEqual(conical_partition.derivative_vertex_angles, expected_derivative_angles);
        testCase.assertEqual(conical_partition.derivative_vertices, expected_derivative_vertices);
        testCase.assertEqual(conical_partition.state_vertex_angles, expected_state_angles);
        testCase.assertEqual(conical_partition.state_vertices, expected_state_vertices);
    end % End of function.

    % ╭──────────────────────────────────────────────╮
    % │             normalizeIndex tests             │
    % ╰──────────────────────────────────────────────╯
    function test_normalizeIndex_2D_single(testCase)
        % ⋘────────── Setup ───────────⋙
        A = eye(2);
        n_slices = 40;
        conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
        
        % ⋘────────── Execute and Verify─────────⋙
        testCase.assertEqual(conical_partition.normalizeIndex(1), 1);
        testCase.assertEqual(conical_partition.normalizeIndex(2), 2);
        testCase.assertEqual(conical_partition.normalizeIndex(n_slices), n_slices);
        testCase.assertEqual(conical_partition.normalizeIndex(n_slices + 1), 1);

    end % End of function.

    % ╭───────────────────────────────────────────────────╮
    % │             getBoundaryVertices tests             │
    % ╰───────────────────────────────────────────────────╯
    function test_getBoundaryVertices_2D(testCase)
      % ⋘────────── Setup ───────────⋙
      A = eye(2);
      n_slices = 30;
      conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
    
      % ⋘────────── Execute and Verify ──────────⋙
      testCase.assertEqual(conical_partition.getBoundaryVertices(1), [1, 2]);
      testCase.assertEqual(conical_partition.getBoundaryVertices(30), [30, 1]);
      testCase.assertEqual(conical_partition.getBoundaryVertices(0), [30, 1]);
      
    end % End function.
    
    % ╭───────────────────────────────────────────────────────────╮
    % │             getConjoiningRegionsIndices tests             │
    % ╰───────────────────────────────────────────────────────────╯
    function test_getConjoiningRegionsIndices_2D(testCase)
      % ⋘────────── Setup ───────────⋙
      A = eye(2);
      n_slices = 30;
      conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
    
      % ⋘────────── Execute and Verify ──────────⋙
      testCase.assertEqual(conical_partition.getConjoiningRegionsIndices(1, 2), 1);
      testCase.assertEqual(conical_partition.getConjoiningRegionsIndices(4, 5), 4);% Check in the middle
      testCase.assertEqual(conical_partition.getConjoiningRegionsIndices(1, 30), 30);% Check at the end, looping to start
      testCase.assertEqual(conical_partition.getConjoiningRegionsIndices(2, 1), 1);% Check with order of given ndxs reversed
    end % End test case function.

    function test_getConjoiningRegionsIndices_emptyWhenNdxsNotAdjacent(testCase)
      % ⋘────────── Setup ───────────⋙
      A = eye(2);
      n_slices = 30;
      conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
    
      % ⋘────────── Execute and Verify ──────────⋙
      testCase.assertEmpty(conical_partition.getConjoiningRegionsIndices(1, 3));
    end % End test case function.
    
    
    % ╭────────────────────────────────────────────────────────────╮
    % │             getStateRegionIndexContainingPoint             │
    % ╰────────────────────────────────────────────────────────────╯
    function test_getStateRegionIndexContainingPoint_2D_atOrigin(testCase)
      % ⋘────────── Setup ───────────⋙
      A = eye(2);
      n_slices = 30;
      conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
      
      % ⋘────────── Execute ─────────⋙
      region_ndxs = conical_partition.getStateRegionIndexContainingPoint([0; 0]);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(region_ndxs, 1:n_slices);
    end % End of function.
    
    
    function test_getStateRegionIndexContainingPoint_2D_notAtOrigin(testCase)
      % ⋘────────── Setup ───────────⋙
      A = magic(2);
      n_slices = 4;
      conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
    
      % Pick xdot (in the derivative space) and x (in the state space) that should be in the cones with index "1". 
      derivative_cone_angle = (2*pi / n_slices);
      for expected_quadrant = 1:n_slices
        % Pick x_dot to to be halfway between the vertices.
        x_dot_angle =  derivative_cone_angle * (expected_quadrant - 1/2);
        xdot = [cos(x_dot_angle); sin(x_dot_angle)];
    
        % Map from xdot to x.
        x    = A \ xdot;
        
        % ⋘────────── Execute ─────────⋙
        state_region_ndxs = conical_partition.getStateRegionIndexContainingPoint(x);
        derivative_region_ndxs = conical_partition.getDerivativeRegionIndexContainingPoint(xdot);
        
        % ⋘────────── Verify ──────────⋙
        testCase.assertEqual(derivative_region_ndxs, expected_quadrant);
        testCase.assertEqual(state_region_ndxs, expected_quadrant);
      end
    
    end % End of function.
    
    % ╭────────────────────────────────────────────────────╮
    % │             getStateRegionIntersections            │
    % ╰────────────────────────────────────────────────────╯
    function test_getStateRegionIntersections_2D_inSingleStateRegion(testCase)
      % ⋘────────── Setup ───────────⋙
      A = eye(2); % Using the identity makes the derivative and state regions equal.
      n_slices = 4;  % Use quadrants.
      conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
      poly = ConvexPolyhedron.fromConvexHull([[1; 1], [0; 0], [1; 1/2]]);
    
      % ⋘────────── Execute ─────────⋙
      [region_ndxs, poly_in_regions] = conical_partition.getStateRegionIntersections(poly);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(region_ndxs, 1);
      testCase.assertNumElements(poly_in_regions, n_slices);
      testCase.assertTrue(poly_in_regions{1} == poly);
      testCase.assertEmpty(poly_in_regions{2});
      testCase.assertEmpty(poly_in_regions{3});
      testCase.assertEmpty(poly_in_regions{4});
    end % End of function.
    
    function test_getStateRegionIntersections_2D_inTwoStateRegion(testCase)
      % ⋘────────── Setup ───────────⋙
      A = eye(2); % Using the identity makes the derivative and state regions equal.
      n_slices = 4;  % Use quadrants.
      conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
    
      % Make "poly" an upward facing wedge, like "\/".
      poly = ConvexPolyhedron.fromConvexHull([[-1; 1], [0; 0], [1; 1]]);
    
      % ⋘────────── Execute ─────────⋙
      [region_ndxs, poly_in_regions] = conical_partition.getStateRegionIntersections(poly);
      
      % ⋘────────── Verify ──────────⋙
      % The intersection of "poly" with the first quadrand should be like "|/".
      expected_poly_1 = ConvexPolyhedron.fromConvexHull([[ 0; 1], [0; 0], [1; 1]]);
      % The intersection of "poly" with the second quadrand should be like "\|".
      expected_poly_2 = ConvexPolyhedron.fromConvexHull([[-1; 1], [0; 0], [0; 1]]);
      testCase.assertEqual(region_ndxs, [1, 2]);
      testCase.assertTrue(poly_in_regions{1} == expected_poly_1);
      testCase.assertTrue(poly_in_regions{2} == expected_poly_2);
    end % End of function.
    
    
    % ! This next test is important as a test of distinct behavior from test_getStateRegionIntersections_2D_inTwoStateRegion because it will cause an entire region to be inside the given polyhedron.
    function test_getStateRegionIntersections_2D_inThreeStateRegion(testCase)
      % ⋘────────── Setup ───────────⋙
      A = eye(2); % Using the identity makes the derivative and state regions equal.
      n_slices = 4;  % Use quadrants.
      conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
    
      % Make "poly" contains all of quadrant 1, along with parts of 2 and 4. 
      poly = ConvexPolyhedron.fromConvexHull([[-1; 1], [0; 0], [1; 1], [1; -0.5]]);
    
      % ⋘────────── Execute ─────────⋙
      [region_ndxs, poly_in_regions] = conical_partition.getStateRegionIntersections(poly);
      
      % ⋘────────── Verify ──────────⋙
      expected_poly_1 = ConvexPolyhedron.fromConvexHull([[ 0; 1], [0; 0], [1; 0], [1; 1]]);
      expected_poly_2 = ConvexPolyhedron.fromConvexHull([[-1; 1], [0; 0], [0; 1]]);
      expected_poly_4 = ConvexPolyhedron.fromConvexHull([[1; 0], [0; 0], [1; -0.5]]);
      testCase.assertEqual(region_ndxs, [1, 2, 4]);
      testCase.assertTrue(poly_in_regions{1} == expected_poly_1, "poly_in_regions{1}");
      testCase.assertTrue(poly_in_regions{2} == expected_poly_2, "poly_in_regions{2}");
      testCase.assertEmpty(poly_in_regions{3});
      testCase.assertTrue(poly_in_regions{4} == expected_poly_4, "poly_in_regions{4}");
    end % End of function.

    % ╭─────────────────────────────────────────────────────────────╮
    % │             getRegionsBetweenVerticesFromAngles             │
    % ╰─────────────────────────────────────────────────────────────╯
    function test_getRegionsBetweenVerticesFromAngles_oneRegion_noWrapping(testCase)
      % ⋘────────── Setup ───────────⋙
      A = eye(2); % Using the identity makes the derivative and state regions equal.
      n_slices = 4;  % Use quadrants.
      conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
      
      % ⋘────────── Execute and Verify ─────────⋙
      testCase.assertEqual(conical_partition.getRegionsBetweenVerticesFromAngles(0, pi/2), 1);
      testCase.assertEqual(conical_partition.getRegionsBetweenVerticesFromAngles(pi/2, pi), 2);
      testCase.assertEqual(conical_partition.getRegionsBetweenVerticesFromAngles(pi, 3*pi/2), 3);
      testCase.assertEqual(conical_partition.getRegionsBetweenVerticesFromAngles(3*pi/2, 2*pi), 4);
    end % End of function.

    function test_getRegionsBetweenVerticesFromAngles_twoRegions_noWrapping(testCase)
      % ⋘────────── Setup ───────────⋙
      A = eye(2); % Using the identity makes the derivative and state regions equal.
      n_slices = 4;  % Use quadrants.
      conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
      
      % ⋘────────── Execute and Verify ─────────⋙
      testCase.assertEqual(conical_partition.getRegionsBetweenVerticesFromAngles(0, pi), [1, 2]);
      testCase.assertEqual(conical_partition.getRegionsBetweenVerticesFromAngles(pi, 2*pi), [3, 4]);
    end % End of function.

    function test_getRegionsBetweenVerticesFromAngles_wrappedPast2pi(testCase)
      % ⋘────────── Setup ───────────⋙
      A = eye(2); % Using the identity makes the derivative and state regions equal.
      n_slices = 4;  % Use quadrants.
      conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
      
      % ⋘────────── Execute and Verify ─────────⋙
      testCase.assertEqual(conical_partition.getRegionsBetweenVerticesFromAngles(3*pi/2, pi/2), [4, 1]);
    end % End of function.

    % ╭──────────────────────────────────────────────╮
    % │ ╭──────────────────────────────────────────╮ │
    % │ │             Static Functions             │ │
    % │ ╰──────────────────────────────────────────╯ │
    % ╰──────────────────────────────────────────────╯
    
    % ╭─────────────────────────────────────────────────╮
    % │             validateAngles function             │
    % ╰─────────────────────────────────────────────────╯
    function test_validateAngles_OK(testCase)
      % ⋘────────── OK ───────────⋙
      testCase.assertWarningFree(@() ConicalPartition.validateAngles([0, 1, 2, 3, 5]));
      testCase.assertWarningFree(@() ConicalPartition.validateAngles([0, 1, 2, 3, 5], [pi-0.1, 1, 2, 3, 5]));
    end

    function test_validateAngles_errorIfDerivativeErrorNondecreasing(testCase)
      errID = "ConicalPartition:DerivativeAnglesNonincreasing";
      testCase.assertError(@() ConicalPartition.validateAngles([0, 1, 1, 2, 3]), errID);
      testCase.assertError(@() ConicalPartition.validateAngles([0, 1, 0.5, 2, 3]), errID);
    end % End of function.

    function test_validateAngles_errorIfAnglesTooBig(testCase)
      errID = "ConicalPartition:DerivativeAngleStepTooBig";
      testCase.assertError(@() ConicalPartition.validateAngles([0, pi+0.1]), errID);
    end % End of function.

    function test_validateAngles_errorIfDerivativeAnglesOutOfRange(testCase)
      errID = "ConicalPartition:DerivativeAnglesOutOfRange";
      testCase.assertError(@() ConicalPartition.validateAngles([0, 1, 2, 3, 7]), errID);
      testCase.assertError(@() ConicalPartition.validateAngles([-2, 0, 1, 2, 3, 4, 5]), errID);
      testCase.assertError(@() ConicalPartition.validateAngles([-2, 0, 1, 2, 3, 8]), errID);
    end % End of function.

    
    % ! We don't need to test for state angles > 180 degrees because this is implied the linear relationship with the derivative angles and the fact that those angles are less than 180 degres.
    % function test_validateAngles_errorIfStateAnglesOutOfRange(testCase)
    % end % End of function.

  end % End methods

end

