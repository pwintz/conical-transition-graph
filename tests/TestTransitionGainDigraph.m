classdef TestTransitionGainDigraph < matlab.unittest.TestCase
% ! You can run these tests via 
% ` runtests TestTransitionGainDigraph

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
      testCase.assertEqual(tg_graph.numVertexNodes(1), cp.n_vertices);
      testCase.assertEqual(tg_graph.numConeNodes(1), cp.n_cones);
      testCase.assertEqual(tg_graph.numNodes(), cp.n_cones + cp.n_vertices);
      testCase.assertEqual(tg_graph.numEdges(), 0);


      % Check that "hasEdgesFromConeVertex" will return false if there is not an edge.
      mode_ndx = 1;
      testCase.assertFalse(tg_graph.hasEdgeFromVertexToCone(  mode_ndx, 2, mode_ndx, 2));
      testCase.assertFalse(tg_graph.hasEdgeFromConeToCone(    mode_ndx, 2, mode_ndx, 2));
      testCase.assertFalse(tg_graph.hasEdgeFromVertexToVertex(mode_ndx, 2, mode_ndx, 2));
      testCase.assertFalse(tg_graph.hasEdgeFromConeToVertex(  mode_ndx, 2, mode_ndx, 2));
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

      % ╭──────────────────────────────────────────────╮
      % │  ╭────────────────────────────────────────╮  │
      % │  │             Multiple Modes             │  │
      % │  ╰────────────────────────────────────────╯  │
      % ╰──────────────────────────────────────────────╯
      function test_multiple_modes(testCase)
        % ⋘────────── Setup ───────────⋙
        conical_partitions = {
          ConicalPartition.fromNumSlices(4);
          ConicalPartition.fromNumSlices(8);
        }
        transition_graph = TransitionGainDigraph(conical_partitions{:});

        mode_1 = 1;
        cone_1 = 2;
        mode_2 = 2;
        cone_2 = 5;
        min_gain = 0.2;
        max_gain = 2.2;

        % ⋘────────── Execute ─────────⋙
        transition_graph.addEdgeFromConeToCone(mode_1, cone_1, mode_2, cone_2, min_gain, max_gain);
        
        % ⋘────────── Verify ──────────⋙
        [edge_row, start_node, end_node] = transition_graph.getEdgeRow(1);
        testCase.assertEqual(start_node.ModeIndex, mode_1);
        testCase.assertEqual(end_node.ModeIndex,   mode_2);
        testCase.assertEqual(start_node.ConeIndex, cone_1);
        testCase.assertEqual(end_node.ConeIndex,   cone_2);
        testCase.assertEqual(start_node.VertexIndex, 0);
        testCase.assertEqual(end_node.VertexIndex,   0);
        testCase.assertEqual(edge_row.MinGain, min_gain);
        testCase.assertEqual(edge_row.MaxGain, max_gain);
        
        % ⋘────────── Setup ───────────⋙
        mode_1 = 2;
        vertex_1 = 2;
        mode_2 = 1;
        cone_2 = 5;
        min_gain = 1.2;
        max_gain = 22.2;

        % ⋘────────── Execute ─────────⋙
        transition_graph.addEdgeFromVertexToVertex(mode_1, vertex_1, mode_2, cone_2, min_gain, max_gain);
        
        % ⋘────────── Verify ──────────⋙
        [edge_row, start_node, end_node] = transition_graph.getEdgeRow(2);
        testCase.assertEqual(start_node.ModeIndex, mode_1);
        testCase.assertEqual(end_node.ModeIndex,   mode_2);
        testCase.assertEqual(start_node.VertexIndex, vertex_1);
        testCase.assertEqual(end_node.VertexIndex,   cone_2);
        testCase.assertEqual(start_node.ConeIndex, 0);
        testCase.assertEqual(end_node.ConeIndex,   0);
        testCase.assertEqual(edge_row.MinGain, min_gain);
        testCase.assertEqual(edge_row.MaxGain, max_gain);
      end % End of function.


      % ╭────────────────────────────────────────────────╮
      % │  ╭──────────────────────────────────────────╮  │
      % │  │             Integration Test             │  │
      % │  ╰──────────────────────────────────────────╯  │
      % ╰────────────────────────────────────────────────╯
      function test_multiple_modes_from_flows_and_jumps(testCase)
        % ⋘────────── Setup ───────────⋙
        flow_map_matrix_1 = [
          2, -3.4;
          2, 1;
        ];
        flow_map_matrix_2 = [
          -1, 3.4;
          -3.4, -2;
        ];

        % Flow set angles
        flow_set_start_angle_1 = 0;
        flow_set_end_angle_1   = pi/2;
        flow_set_start_angle_2 = 0;
        flow_set_end_angle_2   = pi;
        
        flow_set_angles_1 = generateFlowBasedConeAngles2D(...
          flowMapMatrix=flow_map_matrix_1, ...
          flowSetStartAngle=flow_set_start_angle_1, ...
          flowSetEndAngle=flow_set_end_angle_1);
        
        flow_set_angles_2 = generateFlowBasedConeAngles2D(...
          flowMapMatrix=flow_map_matrix_2, ...
          flowSetStartAngle=flow_set_start_angle_2, ...
          flowSetEndAngle=flow_set_end_angle_2);
        
        % Discrete jump maps.
        jump_map_matrix_1_mapsto_1 = [-1.6, 3; -0.2, -0.1]; 
        jump_map_matrix_1_mapsto_2 = -magic(2); 
        jump_map_matrix_2_mapsto_1 = magic(2); 
        jump_map_matrix_2_mapsto_2 = -eye(2) + 0.1; 

        jump_set_1_mapsto_1_start_angle  = pi - 0.1;  jump_set_1_mapsto_1_end_angle    = pi + 0.3;
        jump_set_1_mapsto_2_start_angle  = 5*pi/6;      jump_set_1_mapsto_2_end_angle    = pi;
        jump_set_2_mapsto_1_start_angle  = pi - 0.1;  jump_set_2_mapsto_1_end_angle    = pi;
        jump_set_2_mapsto_2_start_angle  = 0;         jump_set_2_mapsto_2_end_angle    = pi/6;

        jump_set_image_angles_1_mapsto_1 = mapAnglesByLinearMap(jump_map_matrix_1_mapsto_1, [jump_set_1_mapsto_1_start_angle, jump_set_1_mapsto_1_end_angle]);
        jump_set_image_angles_1_mapsto_2 = mapAnglesByLinearMap(jump_map_matrix_1_mapsto_2, [jump_set_1_mapsto_2_start_angle, jump_set_1_mapsto_2_end_angle]);
        jump_set_image_angles_2_mapsto_1 = mapAnglesByLinearMap(jump_map_matrix_2_mapsto_1, [jump_set_2_mapsto_1_start_angle, jump_set_2_mapsto_1_end_angle]);
        jump_set_image_angles_2_mapsto_2 = mapAnglesByLinearMap(jump_map_matrix_2_mapsto_2, [jump_set_2_mapsto_2_start_angle, jump_set_2_mapsto_2_end_angle]);
        
        mode_1_angles = [ ... 
          ... % Add angles for flows.
          flow_set_angles_1, ... 
          ... % Add angles for jumps from this mode.
          jump_set_1_mapsto_1_start_angle, ...
          jump_set_1_mapsto_1_end_angle, ...
          jump_set_1_mapsto_2_start_angle, ...
          jump_set_1_mapsto_2_end_angle, ...
          ... % Add angles for jumps into this mode.
          jump_set_image_angles_1_mapsto_1, ... 
          jump_set_image_angles_2_mapsto_1 ... 
        ];
        
        mode_2_angles = [ ... 
          ... % Add angles for flows.
          flow_set_angles_2, ... 
          ... % Add angles for jumps from this mode.
          jump_set_2_mapsto_1_start_angle, ...
          jump_set_2_mapsto_1_end_angle, ...
          jump_set_2_mapsto_2_start_angle, ...
          jump_set_2_mapsto_2_end_angle, ...
          ... % Add angles for jumps into this mode.  
          jump_set_image_angles_1_mapsto_2, ... 
          jump_set_image_angles_2_mapsto_2 ... 
        ];

        conical_partition_1 = ConicalPartition.fromAngles2D(mode_1_angles);
        conical_partition_2 = ConicalPartition.fromAngles2D(mode_2_angles);

        % ⋘────────── Execute ─────────⋙
        
        
        % ⋘────────── Verify ──────────⋙
        testCase.verifyFail("Test case needs to be implemented.");
      end % End of function.

  end % End of test methods


end
