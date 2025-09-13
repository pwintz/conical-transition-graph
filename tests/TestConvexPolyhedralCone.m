classdef TestConvexPolyhedralCone < matlab.unittest.TestCase
  % ! You can run these tests with 
  % ` runtests TestConvexPolyhedron

  % ╭───────────────────────────────────╮
  % │ ╭───────────────────────────────╮ │
  % │ │             Tests             │ │
  % │ ╰───────────────────────────────╯ │
  % ╰───────────────────────────────────╯
  methods (Test)


  % ╭───────────────────────────────────╮
  % │             From Rays             │
  % ╰───────────────────────────────────╯
  function test_fromRays_2rays_2D(testCase)
    % ⋘────────── Setup ───────────⋙
    
    % ⋘────────── Execute ─────────⋙
    cone = ConvexPolyhedralCone.fromRays([[1; 1], [1; 0]]);
    
    
    % ⋘────────── Verify ──────────⋙
    testCase.assertEqual(cone.n_rays, 2);
    
  end % End of function.

  function test_fromRays_2rays_3D(testCase)
    % ⋘────────── Setup ───────────⋙
    
    % ⋘────────── Execute ─────────⋙
    % Create a 2D cone in 3D space.
    cone = ConvexPolyhedralCone.fromRays([[1; 0; 0], [0; 1; 0]]);
    
    % ⋘────────── Verify ──────────⋙
    testCase.assertEqual(cone.n_rays, 2);
  end % End of function.

  function test_fromRays_3rays_3D(testCase)
    % ⋘────────── Setup ───────────⋙
    
    % ⋘────────── Execute ─────────⋙
    % Create a 2D cone in 3D space.
    cone = ConvexPolyhedralCone.fromRays([[1; 0; 0], [0; 1; 0], [0; 0; 1]]);
    
    % ⋘────────── Verify ──────────⋙
    testCase.assertEqual(cone.n_rays, 3);
  end % End of function.

  function test_fromRays_3rays_3D_not_orthogonal(testCase)
    % ⋘────────── Setup ───────────⋙
    
    % ⋘────────── Execute ─────────⋙
    % Create a 2D cone in 3D space.
    cone = ConvexPolyhedralCone.fromRays([[1; 2; 1], [2; 1; 1], [1; 1; 2]]);
    
    % ⋘────────── Verify ──────────⋙
    testCase.assertEqual(cone.n_rays, 3);
  end % End of function.

    % ╭───────────────────────────────────╮
    % │             isempty()             │
    % ╰───────────────────────────────────╯

    

    % ╭──────────────────────────────────────╮
    % │             Intersection             │
    % ╰──────────────────────────────────────╯
    function test_intersection_cone_cap_cone(testCase)
      % ⋘────────── Setup ───────────⋙
      cone_1 = ConvexPolyhedralCone.fromRays([[1; 0], [1; 2]]);
      cone_2 = ConvexPolyhedralCone.fromRays([[0; 1], [2; 1]]);
      
      % ⋘────────── Execute ─────────⋙
      % Expected rays = [[1; 2], [2; 1]] (normalized)
      p_cap_p = intersection(cone_1, cone_2);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertInstanceOf(p_cap_p, "Polytope");
      testCase.assertTrue(p_cap_p.contains([0; 0]));% Contains origin.
      testCase.assertTrue(p_cap_p.contains([1; 1]));
      testCase.assertTrue(p_cap_p.contains(1e4*[1; 1]));
      testCase.assertFalse(p_cap_p.contains([1; 0]));
      testCase.assertFalse(p_cap_p.contains([0; 1]));
      testCase.assertFalse(p_cap_p.contains([-1; -1]));
    end % End of function.

    function test_intersection_cone_cap_cone_3D(testCase)
      % ⋘────────── Setup ───────────⋙
      cone_1 = ConvexPolyhedralCone.fromRays([[1; 0; 0], [1; 2; 0]]);
      cone_2 = ConvexPolyhedralCone.fromRays([[0; 1; 0], [2; 1; 0]]);
      
      % ⋘────────── Execute ─────────⋙
      % Expected rays = [[1; 2], [2; 1]] (normalized)
      p_cap_p = intersection(cone_1, cone_2);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertInstanceOf(p_cap_p, "Polytope");
      testCase.assertTrue(p_cap_p.contains([0; 0]));% Contains origin.
      testCase.assertTrue(p_cap_p.contains([1; 1]));
      testCase.assertTrue(p_cap_p.contains(1e4*[1; 1]));
      testCase.assertFalse(p_cap_p.contains([1; 0]));
      testCase.assertFalse(p_cap_p.contains([0; 1]));
      testCase.assertFalse(p_cap_p.contains([-1; -1]));
    end % End of function.

    function test_intersection_cone_cap_polytope(testCase)
      % ⋘────────── Setup ───────────⋙
      cone_1 = ConvexPolyhedralCone.fromRays([[1; 0], [1; 2]]);
      cone_2 = Polytope.fromConvexHull([[0; 1], [1; 0], [1; 1], [0; 0]]);
      
      % ⋘────────── Execute ─────────⋙
      % Expected rays = [[1; 2], [2; 1]] (normalized)
      p_cap_p = intersection(cone_1, cone_2);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertInstanceOf(p_cap_p, "Polytope");
      testCase.assertTrue(p_cap_p.contains([0; 0]), "Contains origin");
      testCase.assertTrue(p_cap_p.contains([1; 1]), "Contains [1; 1]");
      testCase.assertFalse(p_cap_p.contains(1e4*[1; 1]), "Does not contain a large vector");
      testCase.assertTrue(p_cap_p.contains([1; 0]));
      testCase.assertFalse(p_cap_p.contains([0; 1]));
      testCase.assertFalse(p_cap_p.contains([-1; -1]));
    end % End of function.
    
  end % End of test methods

end
