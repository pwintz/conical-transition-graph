classdef TestPolytope < matlab.unittest.TestCase
  % ! You can run these tests with 
  % ` runtests TestPolytope

  % ╭───────────────────────────────────╮
  % │ ╭───────────────────────────────╮ │
  % │ │             Tests             │ │
  % │ ╰───────────────────────────────╯ │
  % ╰───────────────────────────────────╯
  methods (Test)

    % ╭──────────────────────────────────────────────────╮
    % │  ╭────────────────────────────────────────────╮  │
    % │  │             Inersection Method             │  │
    % │  ╰────────────────────────────────────────────╯  │
    % ╰──────────────────────────────────────────────────╯
    function test_intersection_polytopes_cap_polytope(testCase)
      % ⋘────────── Setup ───────────⋙
      p1 = Polytope.getBox(2); % [-1, 1] x [-1, 1] box in 2D
      p2 = Polytope.fromTriangle([[1; 0], [0; 0], [0; 1]]);
      
      % ⋘────────── Execute ─────────⋙
      % Expected vertices = [[1; 0], [0; 0], [0; 1]]
      p_cap_p = intersection(p1, p2);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertInstanceOf(p_cap_p, "Polytope");
      testCase.assertTrue(p_cap_p.containsPoints([1/3; 1/3]));
      testCase.assertFalse(p_cap_p.containsPoints([1; 1]));
      testCase.assertFalse(p_cap_p.containsPoints([-1; 1]));
      testCase.assertFalse(p_cap_p.containsPoints([1; -1]));
    end % End of function.
  end % End of test methods

end
