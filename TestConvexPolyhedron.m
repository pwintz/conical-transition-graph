classdef TestConvexPolyhedron < matlab.unittest.TestCase
  % ! You can run these tests to be run automically via ConvexPolyhedron.test() by using the
  % ! "runTestsStaticFunction" snippet to add a static method to ConvexPolyhedron.m

  % methods (TestMethodSetup)
  %	 function setup(testCase)
  %		 % Setup code
  %	 end
  % end
  %
  % methods (TestMethodTeardown)
  %	 function teardown(testCase)
  %		 % Teardown code
  %	 end
  % end

  methods(Static)
    function runTests(varargin) % Define convenience functions for running tests.

      % ⋘──────── Build a list of test strings ────────⋙
      test_class_nam = "TestConvexPolyhedron";
      if isempty(varargin) 
        % If no arguments given, run all of the tests.
        test_strings = test_class_name;
      else
        % If an argument is given, run construct a list of test functions to run.
        test_strings = cellfun(@(function_name) test_class_nam + "/" + function_name, varargin);
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
    % ╭───────────────────────────────────────────────────╮
    % │             Equality "==" overloading             │
    % ╰───────────────────────────────────────────────────╯
    function test_eq_basicEquality(testCase)
      % ⋘────────── Setup ───────────⋙
      cp1_equal    = ConvexPolyhedron.fromConvexHull([[1; 1], [0; 0], [-1; 2]]);
      cp2_equal    = ConvexPolyhedron.fromConvexHull([[1; 1], [0; 0], [-1; 2]]);
      cp3_notEqual = ConvexPolyhedron.fromConvexHull([[1; 1], [0; 3], [-1; 2]]);
      
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
      cp1_equal    = ConvexPolyhedron.fromConvexHull([[1; 1], [0; 0], [-1; 2]]);
      cp2_equal    = ConvexPolyhedron.fromConvexHull([[0; 0], [1; 1], [-1; 2]]);
      cp3_notEqual = ConvexPolyhedron.fromConvexHull([[1; 1], [0; 3], [-1; 2]]);
      
      % ⋘────────── Execute and  Verify ──────────⋙
      % Test operators correctly produce "true".
      testCase.assertTrue( cp1_equal == cp1_equal,    "Same object with ""=="" operator");
      testCase.assertTrue( cp1_equal == cp2_equal,    "Equal with ""=="" operator");
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
      p = ConvexPolyhedron.fromConvexHull([[10; 0], [ 0; 0], [0; 10], [10; 10]]);

      % ⋘────────── Execute and Verify ──────────⋙
      testCase.verifyTrue(p.contains([1; 1]),   "Interior point");
      testCase.verifyTrue(p.contains([10; 0]),  "Vertex");
      testCase.verifyFalse(p.contains([-1; 0]), "Exterior point");
      
      % ! Testing the inclusion of boundary points may be unreilable due to numerical issues, but it should work in this case. 
      testCase.verifyTrue(p.contains([1; 0]), "Boundary point (not vertex)");
    end % End of function.


    function test_contains_multiplePoints(testCase)
      % ⋘────────── Setup ───────────⋙
      p = ConvexPolyhedron.fromConvexHull([[10; 0], [ 0; 0], [0; 10], [10; 10]]);

      % ⋘────────── Execute and Verify ──────────⋙
      testCase.verifyEqual(p.contains([[1; 1], [10; 0], [-1; 0]]), logical([1 1 0]));
    end % End of function.
  end % End of test methods

end
