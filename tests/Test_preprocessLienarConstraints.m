classdef Test_preprocessLienarConstraints < matlab.unittest.TestCase
% ` runtests Test_preprocessLienarConstraints
  methods (Test)
    function test_problematic_case(testCase)
      % ⋘────────── Setup ───────────⋙
      A_ineq=[
              -0.17,      -0.056,         0.98;
              0.17,       -0.056,       -0.98;
              -1.4e-16,   -1,       0;
              1.2e-16,     1,       0;
              0.15,        0.47,       0.87;
              0.46, -0.39, -0.8;
            ]
      b_ineq=[-2.8e-17;
              2.4e-17;
              2.2e-18;
              -5.6e-17;
              1.8e-17;
              5.6e-17]

      [A_ineq, b_ineq, A_eq, b_eq] = preprocessLinearConstraints(A_ineq, b_ineq)

      Ab_ineq = [A_ineq, b_ineq];
      Ab_eq   = [A_eq, b_eq];
      
      % ⋘────────── Execute ─────────⋙
      
      
      % ⋘────────── Verify ──────────⋙
    end % End of function.
  end % End of test methods

end % End of test methods block
