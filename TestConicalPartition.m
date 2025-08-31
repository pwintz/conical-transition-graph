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

  methods(Static)
    function runTests(varargin) % Define convenience functions for running tests.
      % ⋘──────── Build a list of test strings ────────⋙
      test_class_name =  "TestConicalPartition";
      if isempty(varargin) 
        % If no arguments given, run all of the tests.
        test_strings = test_class_name;
      else
        % If an argument is given, run construct a list of test functions to run.
        test_strings = cellfun(@(function_name) test_class_name + "/" + function_name, varargin);
      end
  
      % ⋘──────── Run tests ────────⋙
      results = runtests(test_strings);
  
      % ⋘──────── Print results ────────⋙
      fprintf(...
        "%d Passed, %d Failed, %d Incomplete.\n", ...
        sum([results.Passed]), sum([results.Failed]), sum([results.Incomplete])...
      );
    end % End of function
  end % End static methods block

  % ╭───────────────────────────────────╮
  % │ ╭───────────────────────────────╮ │
  % │ │             Tests             │ │
  % │ ╰───────────────────────────────╯ │
  % ╰───────────────────────────────────╯
  methods (Test)

    % ╭─────────────────────────────────────╮
    % │             Constructor             │
    % ╰─────────────────────────────────────╯
    function test_Constructor_givenAnglesSortedAndUnique(testCase)
      % ⋘────────── Setup ───────────⋙
      angles = 0:0.1:2*pi;
      
      % ⋘────────── Execute ─────────⋙
      conical_partition = ConicalPartition(angles);
      
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(conical_partition.vertex_angles, angles);
      testCase.assertEqual(conical_partition.n_vertices,  numel(angles) + 1); % "+1" for the origin.
      testCase.assertEqual(conical_partition.n_cones,  numel(angles)); 
      
      % testCase.verifyFail("Test case needs to be implemented.");
    end % End of function.
    
    function test_Constructor_errorIfAnglesLargerThanPi(testCase)
      % ⋘────────── Setup ───────────⋙
      angles = [0, pi+0.1, 2*pi - 0.1];

      % ⋘────────── Execute and Verify ─────────⋙
      errID = "ConicalPartition:AngleStepTooBig";
      testCase.assertError(@() ConicalPartition(angles), errID);
    end % End of function.

%     % ╭─────────────────────────────────────────────────────────────╮
%     % │             fromNumberOfUniformDerivativeSlices             │
%     % ╰─────────────────────────────────────────────────────────────╯
%     function test_fromNumberOfUniformDerivativeSlices(testCase)
%         % ⋘────────── Setup ───────────⋙
%         A = [0, -1; 1, 0]; % <- rotation by pi/2
%         n_slices = 4;
%         conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
% 
%         % ⋘────────── Execute and Verify ──────────⋙
%         expected_derivative_vertex_angles = [0, pi/2, pi, 3*pi/2];
%         expected_state_vertex_angles      = mod([0, pi/2, pi, 3*pi/2] - pi/2, 2*pi);
%         testCase.assertEqual(conical_partition.n_vertices, n_slices);
% 
%         derivative_angle_error = pwintz.math.angleDistance(conical_partition.derivative_vertex_angles, expected_derivative_vertex_angles);
%         state_angle_error = pwintz.math.angleDistance(conical_partition.state_vertex_angles, expected_state_vertex_angles);
%         testCase.assertLessThan(derivative_angle_error, 1e-6);
%         testCase.assertLessThan(state_angle_error, 1e-6);
%         
%         testCase.assertEqual(conical_partition.derivative_vertices,     [[1; 0], [0; 1], [-1; 0], [0; -1]], "AbsTol", 1e-6);
%         testCase.assertEqual(conical_partition.state_vertices, [[0; -1], [1; 0], [0; 1], [-1; 0]],          "AbsTol", 1e-6);
%     end % End test 

    % ╭────────────────────────────────────╮
    % │             Constructor             │
    % ╰────────────────────────────────────╯
    function test_Constructor_withoutExtraStateVertexAngles(testCase)
        % ⋘────────── Setup ───────────⋙
        vertex_angles = [0, 2, 4];
        
        % ⋘────────── Execute ─────────⋙
        conical_partition = ConicalPartition(vertex_angles);
        
        % ⋘────────── Verify ──────────⋙
        expected_vertices = [[0; 0], [cos(0); sin(0)], [cos(2); sin(2)], [cos(4); sin(4)]];
        testCase.assertEqual(conical_partition.n_vertices, numel(vertex_angles) + 1);
        testCase.assertEqual(conical_partition.vertices, expected_vertices);
    end % End of function.
    
%     function test_Constructor_withExtraStateVertexAngles_withATheIdentity(testCase)
%         % ⋘────────── Setup ───────────⋙
%         A = eye(2);
%         derivative_vertex_angles  = [0, 2*pi/3, 4*pi/3];
%         extra_state_vertex_angles = [pi, pi/4]; % <- Doesn't need to be sorted.
% 
%         % ⋘────────── Execute ─────────⋙
%         conical_partition = ConicalPartition(A, derivative_vertex_angles, extra_state_vertex_angles);
%         
%         % ⋘────────── Verify ──────────⋙
%         expected_state_angles = [0, pi/4, 2*pi/3, pi, 4*pi/3];
%         expected_state_vertices = pwintz.math.angle2UnitVector(expected_state_angles);
%         % Because we are using the identity matrix for A, the state and derivative vertices are the same.
%         expected_derivative_angles = expected_state_angles;
%         expected_derivative_vertices = expected_state_vertices;
% 
%         testCase.assertEqual(conical_partition.n_vertices, 5);
%         testCase.assertEqual(conical_partition.state_vertex_angles, expected_state_angles, "AbsTol", 1e-6);
%         testCase.assertEqual(conical_partition.derivative_vertex_angles, expected_derivative_angles, "AbsTol", 1e-6);
%         testCase.assertEqual(conical_partition.derivative_vertices, expected_derivative_vertices, "AbsTol", 1e-6);
%         testCase.assertEqual(conical_partition.state_vertices, expected_state_vertices, "AbsTol", 1e-6);
%     end % End of function.

%     function test_Constructor_withExtraStateVertexAngles(testCase)
%         % ⋘────────── Setup ───────────⋙
%         A = magic(2);
%         derivative_vertex_angles  = [0, 5, 2, 4, -2];   % <- Not sorted nor in [0, 2*pi).
%         extra_state_vertex_angles = [pi, pi/4, -1, 12]; % <- Not sorted nor in [0, 2*pi).
%         extra_state_vertices      = pwintz.math.angle2UnitVector(extra_state_vertex_angles);
%         extra_derivative_vertices = A * extra_state_vertices;
%         extra_derivative_angles   = mod(pwintz.math.atan2(extra_derivative_vertices), 2*pi);
% 
%         % ⋘────────── Execute ─────────⋙
%         % expected_derivative_vertices =     [[cos(0); sin(0)], [cos(2); sin(2)], [cos(4); sin(4)]];
%         % expected_state_vertices      = A \ [[cos(0); sin(0)], [cos(2); sin(2)], [cos(4); sin(4)]];
%         
%         conical_partition = ConicalPartition(A, derivative_vertex_angles, extra_state_vertex_angles);
%         
%         % ⋘────────── Verify ──────────⋙
%         expected_derivative_angles = sort(mod([derivative_vertex_angles, extra_derivative_angles], 2*pi));
%         expected_derivative_vertices = pwintz.math.angle2UnitVector(expected_derivative_angles);
%         expected_state_vertices      = A \ expected_derivative_vertices;
%         expected_state_angles   = mod(pwintz.math.atan2(expected_state_vertices), 2*pi);
% 
%         testCase.assertEqual(conical_partition.derivative_vertex_angles, expected_derivative_angles);
%         testCase.assertEqual(conical_partition.derivative_vertices, expected_derivative_vertices);
%         testCase.assertEqual(conical_partition.state_vertex_angles, expected_state_angles);
%         testCase.assertEqual(conical_partition.state_vertices, expected_state_vertices);
%     end % End of function.

    % ╭───────────────────────────────────╮
    % │             getVertex             │
    % ╰───────────────────────────────────╯
    function test_getVertex_canGetAllVertices(testCase)
      % ⋘────────── Setup ───────────⋙
      partition_angles = pwintz.arrays.range(0, 2*pi, "n_values", 20, "includeEnd", false);
      conical_partition = ConicalPartition(partition_angles);
      
      % ⋘────────── Execute and Verify ─────────⋙
      testCase.assertEqual(1:conical_partition.n_vertices, conical_partition.vertex_indices);
      for v_ndx = conical_partition.vertex_indices
        v = conical_partition.getVertex(v_ndx);
        testCase.assertNotEmpty(v);
        testCase.assertSize(v, [2, 1]);
      end
    end % End of function.
    
    function test_getVertex_multiple(testCase)
      % ⋘────────── Setup ───────────⋙
      n_slices = 20;
      conical_partition = ConicalPartition.fromNumSlices(n_slices);
      
      % ⋘────────── Execute ─────────⋙
      v_indices = conical_partition.vertex_indices([1, 1, 2]);
      vertices = conical_partition.getVertex(v_indices);
      
      % ⋘──────── Verify ────────⋙
      testCase.assertEqual(vertices(:, 1), vertices(:, 2));
      testCase.assertSize(vertices, [2, 3]);
    end % End of function.
    
    % ╭───────────────────────────────────────╮
    % │             getCone tests             │
    % ╰───────────────────────────────────────╯
    function test_getCone_single2D(testCase)
      % ⋘────────── Setup ───────────⋙
      n_slices = 20;
      conical_partition = ConicalPartition.fromNumSlices(n_slices);
      
      % ⋘────────── Execute and Verify─────────⋙
      for cone_ndx = conical_partition.cone_indices
        cone = conical_partition.getCone(cone_ndx);
        testCase.assertSize(cone, [1, 1]);
        testCase.assertInstanceOf(cone, ?ConvexPolyhedron);
      end
    end % End of function.

    function test_getCone_multiple2D(testCase)
      % ⋘────────── Setup ───────────⋙
      n_slices = 20;
      conical_partition = ConicalPartition.fromNumSlices(n_slices);
      
      % ⋘────────── Execute and Verify─────────⋙
      cone_ndxs = [1, 2, 2];
      cones = conical_partition.getCone(cone_ndxs);
      testCase.assertSize(cones, size(cone_ndxs));
      testCase.assertInstanceOf(cones, ?ConvexPolyhedron);
    end

%     % ╭──────────────────────────────────────────────╮
%     % │             normalizeIndex tests             │
%     % ╰──────────────────────────────────────────────╯
%     function test_normalizeVertexIndex_2D_single(testCase)
%       % ⋘────────── Setup ───────────⋙
%       n_slices = 40;
%       conical_partition = ConicalPartition.fromNumSlices(n_slices);
%       
%       % ⋘────────── Execute and Verify─────────⋙
%       testCase.assertEqual(conical_partition.normalizeVertexIndex(1), 1);
%       testCase.assertEqual(conical_partition.normalizeVertexIndex(2), 2);
%       testCase.assertEqual(conical_partition.normalizeVertexIndex(n_slices), n_slices);
%       testCase.assertEqual(conical_partition.normalizeVertexIndex(n_slices + 1), 1);
%     end % End of function.
% 
%     function test_normalizeConeIndex_2D_single(testCase)
%       % ⋘────────── Setup ───────────⋙
%       n_slices = 40;
%       conical_partition = ConicalPartition.fromNumSlices(n_slices);
%       
%       % ⋘────────── Execute and Verify─────────⋙
%       testCase.assertEqual(conical_partition.normalizeConeIndex(1), 1);
%       testCase.assertEqual(conical_partition.normalizeConeIndex(2), 2);
%       testCase.assertEqual(conical_partition.normalizeConeIndex(n_slices), n_slices);
%       testCase.assertEqual(conical_partition.normalizeConeIndex(n_slices + 1), 1);
%     end % End of function.

    % ╭──────────────────────────────────────────────────────────────────────────────╮
    % │             getVerticesAdjacentToCone, getConesConjoiningTwoVertices tests             │
    % ╰──────────────────────────────────────────────────────────────────────────────╯
    function test_getVerticesAdjacentToCone_2D(testCase)
      % ⋘────────── Setup ───────────⋙
      n_slices = 4;
      conical_partition = ConicalPartition.fromNumSlices(n_slices);

      % ⋘────────── Execute and Verify ──────────⋙
      origin_ndx = conical_partition.origin_index;
      non_origin_vertex_indices = setdiff(conical_partition.vertex_indices, origin_ndx);
      for v1_ndx = non_origin_vertex_indices
        for v2_ndx = non_origin_vertex_indices% (v1_ndx+1):conical_partition.n_vertices
          % Get cone that contains both points.
          if v1_ndx == v2_ndx 
            continue
          end
          % Find the cones that have both v1 and v2 on their boundaries. 
          cone_ndx = conical_partition.getConesConjoiningTwoVertices(v1_ndx, v2_ndx);
          
          % If v1 and v2 are adjacent in the vertex graph, then cone_ndx should be nonempty.
          if conical_partition.areVerticesAdjacent(v1_ndx, v2_ndx)
            % Get the boundary vertices using getVerticesAdjacentToCone().
            v_ndxs_adjacent_to_cone = conical_partition.getVerticesAdjacentToCone(cone_ndx);
            % Check that the boundary vertices are v1, v2, and the origin.
            expected_ndxs = [v1_ndx, v2_ndx, origin_ndx];
            testCase.assertEqual(sort(v_ndxs_adjacent_to_cone), sort(expected_ndxs), ...
              sprintf("v1_ndx = %d, v2_ndx = %d, cone_ndx=%d, edges=%s, faces=%s", v1_ndx, v2_ndx, cone_ndx, mat2str(conical_partition.graph.edges), mat2str(conical_partition.graph.faces))...
            );
          else
            testCase.assertEmpty(cone_ndx);
          end
        end
      end % end for loop

    end % End function.
    
    % ╭────────────────────────────────────────╮
    % │             getConeNormals             │
    % ╰────────────────────────────────────────╯
    function test_getConeNormals_oneOutputArg(testCase)
      % ⋘────────── Setup ───────────⋙
      n_slices = 4;
      conical_partition = ConicalPartition.fromNumSlices(n_slices);
      for cone_ndx = conical_partition.cone_indices
        
        % ⋘────────── Execute ─────────⋙
        cone_normals = conical_partition.getConeNormals(cone_ndx);
        
        % ⋘────────── Verify Boundary Points ──────────⋙
        % Check that the <x, n> = 0 for x in the boundary and n the corresponding normal vector.
        [~, bnd_cone_vertices] = conical_partition.getVerticesAdjacentToCone(cone_ndx);
        bnd_dot_products = dot(bnd_cone_vertices, cone_normals);
        testCase.assertEqual(bnd_dot_products, [0, 0], "AbsTol", 1e-2);
        
        % ⋘────────── Verify Interior Points ──────────⋙
        % Check that <x, n> > 0 for point x in the interior of the cone and for each normal vector n...
        x = bnd_cone_vertices(:, 1) + bnd_cone_vertices(:, 2); % The sum of the boundary points is an interior point because the cone is convex.
        interior_dot_product = x' * cone_normals;
        testCase.assertGreaterThan(interior_dot_product, 0);
      end
    end % End of function.

    function test_getConeNormals_twoOutputArgs(testCase)
      % ⋘────────── Setup ───────────⋙
      n_slices = 8;
      conical_partition = ConicalPartition.fromNumSlices(n_slices);
      for cone_ndx = conical_partition.cone_indices
        
        % ⋘────────── Execute ─────────⋙
        [cone_normals, bnd_vertex_ndxs] = conical_partition.getConeNormals(cone_ndx);
        
        % ⋘────────── Verify Boundary Points ──────────⋙
        % Check that the <x, n> = 0 for x in the boundary and n the corresponding normal vector.
        bnd_vertices = conical_partition.getVertex(bnd_vertex_ndxs);
        bnd_dot_products = dot(bnd_vertices, cone_normals);
        testCase.assertEqual(bnd_dot_products, [0, 0], "AbsTol", 1e-2);
        
        % ⋘────────── Verify Interior Points ──────────⋙
        % Check that <x, n> > 0 for a point x in the interior of the cone and for each normal vector n...
        % The sum of the boundary points is an interior point because the cone is convex.
        x = bnd_vertices(:, 1) + bnd_vertices(:, 2); 
        interior_dot_product = x' * cone_normals;
        testCase.assertGreaterThan(interior_dot_product, 0);
      end
    end % End of function.

    % ╭────────────────────────────────────────────────╮
    % │             getConesContainingPoint            │
    % ╰────────────────────────────────────────────────╯
    function test_getConesContainingPoint_2D_atOrigin(testCase)
      % ⋘────────── Setup ───────────⋙
      angles = 0:0.1:2*pi;
      conical_partition = ConicalPartition(angles);
      
      % ⋘────────── Execute ─────────⋙
      cone_ndxs = conical_partition.getConesContainingPoint([0; 0]);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(cone_ndxs, 1:conical_partition.n_cones);
    end % End of function.
    
    
    function test_getConesContainingPoint_2D_notAtOrigin(testCase)
      % ⋘────────── Setup ───────────⋙
      cone_angle = pi/2;
      partition_angles = cone_angle * (0:3);
      conical_partition = ConicalPartition(partition_angles);
    
      % Pick xdot (in the derivative space) and x (in the state space) that should be in the cones with index "1". 
      for i = 1:(numel(partition_angles))
        % Pick x to to be halfway between the vertices.
        x_angle = partition_angles(i) + 0.5 * cone_angle;
        x = [cos(x_angle); sin(x_angle)];
    
        % ⋘────────── Execute ─────────⋙
        cone_ndxs = conical_partition.getConesContainingPoint(x);
        
        % ⋘────────── Verify ──────────⋙
        testCase.assertEqual(cone_ndxs, i);
      end
    end % End of function.
    
    % ╭─────────────────────────────────────────────╮
    % │             getConeIntersections            │
    % ╰─────────────────────────────────────────────╯
    function test_getConeIntersections_2D_inSingleCone(testCase)
      % ⋘────────── Setup ───────────⋙
      n_slices = 4;
      conical_partition = ConicalPartition.fromNumSlices(n_slices);

      % Create a polyhedron that is wedge-shaped, in the first quadrant.
      polyhedron = ConvexPolyhedron.fromConvexHull([[1; 1], [0; 0], [1; 1/2]]);
    
      % ⋘────────── Execute ─────────⋙
      [intersecting_cone_ndxs, does_cone_intersect, intersections] = conical_partition.getConeIntersections(polyhedron);
      
      % ⋘────────── Verify ──────────⋙
      conical_partition.cones{:}
      expected_cone_ndx = conical_partition.getConesContainingPoint([1; 1])
      conical_partition.getCone(expected_cone_ndx)
      testCase.assertEqual(intersecting_cone_ndxs, expected_cone_ndx);
      testCase.assertNumElements(intersections, n_slices);
      
      testCase.assertTrue(intersections{1} == expected_poly);
      testCase.assertEqual(intersections{1},  expected_poly);
      testCase.assertEmpty(poly_in_cones{2});
      testCase.assertEmpty(poly_in_cones{3});
      testCase.assertEmpty(poly_in_cones{4});
    end % End of function.
    
    function test_getConeIntersections_2D_inTwoCone(testCase)
      % ⋘────────── Setup ───────────⋙
      partition_angles = pwintz.arrays.range(0, 2*pi, "n_values", 4, "includeEnd", false);
      conical_partition = ConicalPartition(partition_angles);
    
      % Make "poly" an upward facing wedge, like "\/".
      poly = ConvexPolyhedron.fromConvexHull([[-1; 1], [0; 0], [1; 1]]);
    
      % ⋘────────── Execute ─────────⋙
      [cone_ndxs, poly_in_cones] = conical_partition.getConeIntersections(poly);
      
      % ⋘────────── Verify ──────────⋙
      % The intersection of "poly" with the first quadrand should be like "|/".
      expected_poly_1 = ConvexPolyhedron.fromConvexHull([[ 0; 1], [0; 0], [1; 1]]);
      % The intersection of "poly" with the second quadrand should be like "\|".
      expected_poly_2 = ConvexPolyhedron.fromConvexHull([[-1; 1], [0; 0], [0; 1]]);
      testCase.assertEqual(cone_ndxs, [1, 2]);
      testCase.assertTrue(poly_in_cones{1} == expected_poly_1);
      testCase.assertTrue(poly_in_cones{2} == expected_poly_2);
    end % End of function.
    
    
    % ! This next test is important as a test of distinct behavior from test_getConeIntersections_2D_inTwoCone because it will cause an entire cone to be inside the given polyhedron.
    function test_getConeIntersections_2D_inThreeCones(testCase)
      % ⋘────────── Setup ───────────⋙
      n_slices = 4;
      conical_partition = ConicalPartition.fromNumSlices(n_slices);

      % Make "poly" contains all of quadrant 1, along with parts of 2 and 4. 
      poly = ConvexPolyhedron.fromConvexHull([[-1; 1], [0; 0], [1; 1], [1; -0.5]]);
    
      % ⋘────────── Execute ─────────⋙
      [cone_ndxs, poly_in_cones] = conical_partition.getConeIntersections(poly);
      
      % ⋘────────── Verify ──────────⋙
      expected_poly_1 = ConvexPolyhedron.fromConvexHull([[ 0; 1], [0; 0], [1; 0], [1; 1]]);
      expected_poly_2 = ConvexPolyhedron.fromConvexHull([[-1; 1], [0; 0], [0; 1]]);
      expected_poly_4 = ConvexPolyhedron.fromConvexHull([[1; 0], [0; 0], [1; -0.5]]);
      testCase.assertEqual(cone_ndxs, [1, 2, 4]);
      testCase.assertTrue(poly_in_cones{1} == expected_poly_1, "poly_in_cones{1}");
      testCase.assertTrue(poly_in_cones{2} == expected_poly_2, "poly_in_cones{2}");
      testCase.assertEmpty(poly_in_cones{3});
      testCase.assertTrue(poly_in_cones{4} == expected_poly_4, "poly_in_cones{4}");
    end % End of function.

    % ╭─────────────────────────────────────────────────────────────╮
    % │             getConesBetweenVerticesConstructor             │
    % ╰─────────────────────────────────────────────────────────────╯
    function test_getConesBetweenVerticesFromAngles_oneCone_noWrapping(testCase)
      % ⋘────────── Setup ───────────⋙
      n_slices = 4;
      conical_partition = ConicalPartition.fromNumSlices(n_slices);
      
      % ⋘────────── Execute and Verify ─────────⋙
      % First quadrant.
      actual_cone_ndx   = conical_partition.getConesIntersectingArc(0, pi/2);
      expected_cone_ndx = conical_partition.getConesContainingPoint([1; 1]);
      testCase.assertEqual(actual_cone_ndx, expected_cone_ndx);

      % Second quadrant.
      actual_cone_ndx   = conical_partition.getConesIntersectingArc(pi/2, pi);
      expected_cone_ndx = conical_partition.getConesContainingPoint([-1; 1]);
      testCase.assertEqual(actual_cone_ndx, expected_cone_ndx);

      % Third quadrant.
      [actual_cone_ndx, arc]   = conical_partition.getConesIntersectingArc(pi, -pi/2);
      expected_cone_ndx = conical_partition.getConesContainingPoint([-1; -1]);
      testCase.assertEqual(actual_cone_ndx, expected_cone_ndx);

      % Fourth quadrant.
      expected_cone_ndx = conical_partition.getConesContainingPoint([1; -1]);
      actual_cone_ndx   = conical_partition.getConesIntersectingArc(3*pi/2, 2*pi);
      testCase.assertEqual(actual_cone_ndx, expected_cone_ndx);

    end % End of function.

    function test_getConesBetweenVerticesFromAngles_twoCones_noWrapping(testCase)
      % ⋘────────── Setup ───────────⋙
      A = eye(2); % Using the identity makes the derivative and state cones equal.
      n_slices = 4;  % Use quadrants.
      conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
      
      % ⋘────────── Execute and Verify ─────────⋙
      testCase.assertEqual(conical_partition.getConesBetweenVerticesFromAngles(0, pi), [1, 2]);
      testCase.assertEqual(conical_partition.getConesBetweenVerticesFromAngles(pi, 2*pi), [3, 4]);
    end % End of function.

    function test_getConesBetweenVerticesFromAngles_wrappedPast2pi(testCase)
      % ⋘────────── Setup ───────────⋙
      A = eye(2); % Using the identity makes the derivative and state cones equal.
      n_slices = 4;  % Use quadrants.
      conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(A, n_slices);
      
      % ⋘────────── Execute and Verify ─────────⋙
      testCase.assertEqual(conical_partition.getConesBetweenVerticesFromAngles(3*pi/2, pi/2), [4, 1]);
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
      errID = "ConicalPartition:AnglesNonincreasing";
      testCase.assertError(@() ConicalPartition.validateAngles([0, 1, 1, 2, 3]), errID);
      testCase.assertError(@() ConicalPartition.validateAngles([0, 1, 0.5, 2, 3]), errID);
    end % End of function.

    function test_validateAngles_errorIfAnglesTooBig(testCase)
      errID = "ConicalPartition:AngleStepTooBig";
      testCase.assertError(@() ConicalPartition.validateAngles([0, pi+0.1]), errID);
    end % End of function.

    function test_validateAngles_errorIfDerivativeAnglesOutOfRange(testCase)
      errID = "ConicalPartition:AnglesOutOfRange";
      testCase.assertError(@() ConicalPartition.validateAngles([0, 1, 2, 3, 7]), errID);
      testCase.assertError(@() ConicalPartition.validateAngles([-2, 0, 1, 2, 3, 4, 5]), errID);
      testCase.assertError(@() ConicalPartition.validateAngles([-2, 0, 1, 2, 3, 8]), errID);
    end % End of function.

    
    % ! We don't need to test for state angles > 180 degrees because this is implied the linear relationship with the derivative angles and the fact that those angles are less than 180 degres.
    % function test_validateAngles_errorIfStateAnglesOutOfRange(testCase)
    % end % End of function.

  end % End methods

end

