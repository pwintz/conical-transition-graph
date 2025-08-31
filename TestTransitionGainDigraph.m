classdef TestTransitionGainDigraph < matlab.unittest.TestCase
% ! You can run these tests to be run automically via TransitionGainDigraph.test() by using the
% ! "runTestsStaticFunction" snippet to add a static method to TransitionGainDigraph.m,
%or by running TestTransitionGainDigraph.runTests()

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
    test_class_name =  "TestTransitionGainDigraph";
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
    function test_Constructor(testCase)
      % ⋘────────── Setup ───────────⋙
      n_slices = 4;
      cp = ConicalPartition.fromNumSlices(n_slices);
      
      % ⋘────────── Execute ─────────⋙
      tg_graph = TransitionGainDigraph(cp);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(tg_graph.numVertexNodes, cp.n_vertices);
      testCase.assertEqual(tg_graph.numConeNodes, cp.n_cones);
      testCase.assertEqual(tg_graph.numNodes, cp.n_cones + cp.n_vertices);
      testCase.assertEqual(tg_graph.numEdges, 0);


      % Check that "hasEdgesFromConeVertex" will return false if there is not an edge.
      testCase.assertFalse(tg_graph.hasEdgeFromVertexToCone(2, 2));
      testCase.assertFalse(tg_graph.hasEdgeFromConeToCone(2, 2));
      testCase.assertFalse(tg_graph.hasEdgeFromVertexToVertex(2, 2));
      testCase.assertFalse(tg_graph.hasEdgeFromConeToVertex(2, 2));
    end % End of function.

    % ╭────────────────────────────────────────────────────────────╮
    % │  ╭──────────────────────────────────────────────────────╮  │
    % │  │             Test Add/Get Edges Functions             │  │
    % │  ╰──────────────────────────────────────────────────────╯  │
    % ╰────────────────────────────────────────────────────────────╯
    function test_edgeFromVertexToCone(testCase)
      % ⋘────────── Setup ───────────⋙
      cp = ConicalPartition.fromNumSlices(4);
      tg_graph = TransitionGainDigraph(cp);
      start_vertex_ndx = 1;
      end_cone_ndx     = 1;
      min_gain = 1.0;
      max_gain = 2.0;
      
      % ⋘────────── Execute ─────────⋙
      tg_graph.addEdgeFromVertexToCone(start_vertex_ndx, end_cone_ndx, min_gain, max_gain);
      has_edge = tg_graph.hasEdgeFromVertexToCone(start_vertex_ndx, end_cone_ndx);
      edge    = tg_graph.getEdgesFromVertexToCone(start_vertex_ndx, end_cone_ndx);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(tg_graph.numEdges, 1);
      testCase.assertTrue(has_edge);
      testCase.assertEqual(edge.MaxGain, max_gain);
      testCase.assertEqual(edge.MinGain, min_gain);
    end % End of function.

    function test_edgeFromVertexToVertex(testCase)
      % ⋘────────── Setup ───────────⋙
      cp = ConicalPartition.fromNumSlices(4);
      tg_graph = TransitionGainDigraph(cp);
      start_vertex_ndx = 1;
      end_vertex_ndx = 1;
      max_gain = 2.0;
      min_gain = 2.0;
      
      % ⋘────────── Execute ─────────⋙
      tg_graph.addEdgeFromVertexToVertex(start_vertex_ndx, end_vertex_ndx, min_gain, max_gain);
      has_edge = tg_graph.hasEdgeFromVertexToVertex(start_vertex_ndx, end_vertex_ndx);
      edge    = tg_graph.getEdgesFromVertexToVertex(start_vertex_ndx, end_vertex_ndx);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(tg_graph.numEdges, 1);
      testCase.assertTrue(has_edge);
      testCase.assertEqual(edge.MaxGain, max_gain);
      testCase.assertEqual(edge.MinGain, min_gain);
    end % End of function.
    
    function test_edgeFromConeToCone(testCase)
      % ⋘────────── Setup ───────────⋙
      cp = ConicalPartition.fromNumSlices(4);
      tg_graph = TransitionGainDigraph(cp);
      start_cone_ndx = 1;
      end_cone_ndx = 1;
      max_gain = 2.0;
      min_gain = 2.0;
      
      % ⋘────────── Execute ─────────⋙
      tg_graph.addEdgeFromConeToCone(start_cone_ndx, end_cone_ndx, min_gain, max_gain);
      has_edge = tg_graph.hasEdgeFromConeToCone(start_cone_ndx, end_cone_ndx);
      edge    = tg_graph.getEdgesFromConeToCone(start_cone_ndx, end_cone_ndx);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(tg_graph.numEdges, 1);
      testCase.assertTrue(has_edge);
      testCase.assertEqual(edge.MaxGain, max_gain);
      testCase.assertEqual(edge.MinGain, min_gain);
    end % End of function.
    
    function test_edgeFromConeToVertex(testCase)
      % ⋘────────── Setup ───────────⋙
      cp = ConicalPartition.fromNumSlices(4);
      tg_graph = TransitionGainDigraph(cp);
      start_cone_ndx = 4;
      end_vertex_ndx = 5;
      max_gain = 2.0;
      min_gain = 2.0;
      
      % ⋘────────── Execute ─────────⋙
      tg_graph.addEdgeFromConeToVertex(start_cone_ndx, end_vertex_ndx, min_gain, max_gain);
      has_edge = tg_graph.hasEdgeFromConeToVertex(start_cone_ndx, end_vertex_ndx);
      edge    = tg_graph.getEdgesFromConeToVertex(start_cone_ndx, end_vertex_ndx);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(tg_graph.numEdges, 1);
      testCase.assertTrue(has_edge);
      testCase.assertEqual(edge.MaxGain, max_gain);
      testCase.assertEqual(edge.MinGain, min_gain);
    end % End of function.

    % ╭───────────────────────────────────────────╮
    % │  ╭─────────────────────────────────────╮  │
    % │  │             Test Cycles             │  │
    % │  ╰─────────────────────────────────────╯  │
    % ╰───────────────────────────────────────────╯
    function test_cycles_singleCycleWithTwoVertices(testCase)
      % ⋘────────── Setup ───────────⋙
      cp = ConicalPartition.fromNumSlices(4);
      tg_graph = TransitionGainDigraph(cp);
      v1_ndx = 1;
      v2_ndx = 2;
      min_gain_1 = 0.5;
      max_gain_1 = 2.0;
      min_gain_2 = 0.3;
      max_gain_2 = 1.2;
      tg_graph.addEdgeFromVertexToVertex(v1_ndx, v2_ndx, min_gain_1, max_gain_1);
      tg_graph.addEdgeFromVertexToVertex(v2_ndx, v1_ndx, min_gain_2, max_gain_2);
      
      % ⋘────────── Execute ─────────⋙
      cycles = tg_graph.getCycles();
      [cycle_min_gains, cycle_max_gains] = tg_graph.getCycleGains()
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertSize(cycles, [1, 1]);
      testCase.assertSize(cycle_min_gains, [1, 1]);
      testCase.assertSize(cycle_max_gains, [1, 1]);
      
      testCase.assertEqual(cycle_min_gains(1), min_gain_1 * min_gain_2, "AbsTol", 1e-6);
      testCase.assertEqual(cycle_max_gains(1), max_gain_1 * max_gain_2, "AbsTol", 1e-6);
      
      
      % testCase.assertEqual(edge.MaxGain, max_gain);
      % testCase.assertEqual(edge.MinGain, min_gain);
    end % End of function.

    % ╭───────────────────────────────────────────────────╮
    % │             getReachableSetFromVertex             │
    % ╰───────────────────────────────────────────────────╯
    % Run via
    % ` runtests TestTransitionGainDigraph/test_getReachableSetFromVertex_empty
    function test_getReachableSetFromVertex_empty(testCase)
      % ⋘────────── Setup ───────────⋙
      cp = ConicalPartition.fromNumSlices(4);
      tg_graph = TransitionGainDigraph(cp);
      
      % ⋘────────── Execute ─────────⋙
      start_vertex_ndx   = 1;
      vertex_reach_table = tg_graph.getReachableSetFromVertex(start_vertex_ndx);
      
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(size(vertex_reach_table, 1), 1);
      
    end % End of function.
    function test_getReachableSetFromVertex_multiple(testCase)
      % ` runtests TestTransitionGainDigraph/test_getReachableSetFromVertex_multiple
        % ⋘────────── Setup ───────────⋙
        cp = ConicalPartition.fromNumSlices(4);
        tg_graph = TransitionGainDigraph(cp);
        tg_graph.addEdgeFromVertexToVertex(1, 2, 0.2, 1.3);
        tg_graph.addEdgeFromVertexToVertex(1, 3, 0.2, 1.3);
        
        % ⋘────────── Execute ─────────⋙
        start_vertex_ndx   = 1;
        vertex_reach_table = tg_graph.getReachableSetFromVertex(start_vertex_ndx);
        
        % ⋘────────── Verify ──────────⋙
        testCase.assertEqual(size(vertex_reach_table, 1), 1);
        
      end % End of function.
  end % End of test methods


end
