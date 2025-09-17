classdef TestFlowTransitionGainDigraph < matlab.unittest.TestCase
% ` runtests TestFlowTransitionGainDigraph

  methods (Test)
    function test_constructor_2D(testCase)
      % ⋘────────── Setup ───────────⋙
      conical_partition = ConicalPartition.fromNumSlices(4);
      flow_set_cone_ndxs = 1:(conical_partition.n_cones / 2);
      flow_map_matrix    = -eye(2); % Points move straight toward origin. 
      
      % ⋘────────── Execute ─────────⋙
      FlowTransitionGainDigraph(conical_partition, flow_set_cone_ndxs, flow_map_matrix);
    end % End of function.

    function test_rotating_flow(testCase)
      % ⋘────────── Setup ───────────⋙
      n_cones = 8;
      conical_partition = ConicalPartition.fromNumSlices(n_cones);
      flow_set_cone_ndxs = 1:n_cones;

      % Rotational vector field
      flow_map_matrix    = pwintz.linear_algebra.rotation2(pi/2);   
      
      % ⋘────────── Execute ─────────⋙
      flow_digraph = FlowTransitionGainDigraph(conical_partition, flow_set_cone_ndxs, flow_map_matrix);
      
      % ⋘────────── Verify ──────────⋙
      % There should be three edges for each cone: 
      % * From the cone to the next vertex
      % * From the previous vertex to the cone
      % * Between the vertices.
      testCase.assertEqual(flow_digraph.numEdges(), 3 * n_cones);

      testCase.assertEqual(size(flow_digraph.getEdgesFromVerticesToVertices(), 1), n_cones);
      testCase.assertEqual(size(flow_digraph.getEdgesFromConesToVertices(), 1), n_cones);
      testCase.assertEqual(size(flow_digraph.getEdgesFromVerticesToCones(), 1), n_cones);

      % testCase.assertNotEmpty(flow_digraph.restricted_reachable_sets_from_unit_sphere_in_cone);
      % testCase.assertNotEmpty(flow_digraph.can_flow_from_vertex_into_cone);
      % testCase.assertNotEmpty(flow_digraph.can_flow_from_cone_to_vertex);

      cones_transition_graph = flow_digraph.contractEdgesToBetweenCones(1:n_cones, 1:n_cones);
      testCase.assertEqual(cones_transition_graph.numEdges(), n_cones^2);

      % TODO: Test cycles
      % cycles = flow_digraph.getCycles();
      % [cycle_min_gains, cycle_max_gains, cycles_nodes, cycles_edges] = flow_digraph.getCycleGains();
      % testCase.assertNumElements(cycles, 1);
      % testCase.assertEqual(cycles{1}, conical_partition.ray_indices);
      %
      % % Since the flow is purely rotational, the gain of any solution is 1.0, but since we are approximating the system, the max and min gains we compute will be over-conservative, meaning:
      % % `
      % % `   cycle_min_gains ≤ 1.0 ≤ cycle_max_gains
      % % `
      % testCase.assertLessThanOrEqual(cycle_min_gains, 1.0);
      % testCase.assertGreaterThanOrEqual(cycle_max_gains, 1.0);

      % testCase.verifyFail("Test case needs to be implemented.");
    end % End of function.

    function test_two_modes(testCase)
      % ⋘────────── Setup ───────────⋙
      conical_partition_1 = ConicalPartition.fromNumSlices(4);
      conical_partition_2 = ConicalPartition.fromNumSlices(8);

      flow_set_cone_ndxs_1 = 1:3;
      flow_set_cone_ndxs_2 = conical_partition_2.cone_indices;
      flow_map_matrix_1 = [-2, -1; 0, -1];
      flow_map_matrix_2 = [-1, 2; -2, -1];

      % ⋘────────── Execute ─────────⋙
      % flow_specs_1 = FlowSpecification(flow_set_cone_ndxs_1, flow_map_matrix_1);
      % flow_specs_2 = FlowSpecification(flow_set_cone_ndxs_2, flow_map_matrix_2);

%       flow_specs = [flow_specs_1, flow_specs_2];
% 
%       conical_partitions = [conical_partition_1, conical_partition_2];

      flow_graph = FlowTransitionGainDigraph(conical_partition_1, flow_set_cone_ndxs_1, flow_map_matrix_1, conical_partition_2, flow_set_cone_ndxs_2, flow_map_matrix_2);
        
%       % ⋘────────── Verify ──────────⋙
%       testCase.assertEqual(jump_graph.jump_set_cone_ndxs{1, 1}, 1:3);
%       testCase.assertEqual(jump_graph.jump_set_image_cone_ndxs{1, 1}, [2, 3, 4]);% In left half plane. 
%       mode_ndx = 1; 
%       testCase.assertTrue(jump_graph.hasEdgeFromConeToCone(mode_ndx, 1, mode_ndx, 4));
%       testCase.assertTrue(jump_graph.hasEdgeFromConeToCone(mode_ndx, 2, mode_ndx, 3));
%       testCase.assertFalse(jump_graph.hasEdgeFromConeToCone(mode_ndx, 1, mode_ndx, 2));
%       testCase.assertFalse(jump_graph.hasEdgeFromConeToCone(mode_ndx, 3, mode_ndx, 4));
%       testCase.assertFalse(jump_graph.hasEdgeFromConeToCone(mode_ndx, 4, mode_ndx, 3));
% 
% 
%       mode_ndx = 2; 
%       x = [1; 2];
%       cone_ndx_1 = conical_partition_2.getConesContainingPoint( x);
%       cone_ndx_2 = conical_partition_2.getConesContainingPoint(-x);
% 
%       testCase.assertTrue(jump_graph.hasEdgeFromConeToCone(mode_ndx, cone_ndx_1, mode_ndx, cone_ndx_2));
%       testCase.assertTrue(jump_graph.hasEdgeFromConeToCone(mode_ndx, cone_ndx_2, mode_ndx, cone_ndx_1));
%       testCase.assertFalse(jump_graph.hasEdgeFromConeToCone(mode_ndx, 1, mode_ndx, 1));
%       testCase.assertFalse(jump_graph.hasEdgeFromConeToCone(mode_ndx, 1, mode_ndx, 2));
% 
%       % Check that getCycles doesn't produce errors.
%       jump_graph.getCycles();
%       % testCase.assertNumElements(cycles, 1); % One cycle. 
% 
%       % testCase.assertEqual(cycles{1}, [2, 3])
      
    end % End of function.

  end % End of test methods

end % End of test methods block
