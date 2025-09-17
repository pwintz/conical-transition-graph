classdef TestFlowReachabilityAnalyzer < matlab.unittest.TestCase
% ` runtests TestFlowReachabilityAnalyzer

  methods (Test)
    function test_constructor_2D(testCase)
      % ⋘────────── Setup ───────────⋙
      conical_partition = ConicalPartition.fromNumSlices(4);
      flow_set_cone_ndxs = 1:(conical_partition.n_cones / 2);
      flow_map_matrix    = -eye(2); % Points move straight toward origin. 
      
      % ⋘────────── Execute ─────────⋙
      FlowReachabilityAnalyzer(conical_partition, flow_set_cone_ndxs, flow_map_matrix);
    end % End of function.

    function test_rotating_flow(testCase)
      % ⋘────────── Setup ───────────⋙
      n_cones = 8;
      conical_partition = ConicalPartition.fromNumSlices(n_cones);
      flow_set_cone_ndxs = 1:n_cones;

      % Rotational vector field
      flow_map_matrix    = pwintz.linear_algebra.rotation2(pi/2);   
      
      % ⋘────────── Execute ─────────⋙
      reach_analyzer = FlowReachabilityAnalyzer(conical_partition, flow_set_cone_ndxs, flow_map_matrix);
      
    
      % ⋘────────── Verify ──────────⋙
      % There should be three edges for each cone: 
      % * From the cone to the next vertex
      % * From the previous vertex to the cone
      % * Between the vertices.
      % testCase.assertEqual(reach_analyzer.numEdges(), 3 * n_cones);

      % Exactly one vertex and one cone is reachable from each ray.
      for ray_ndx = conical_partition.ray_indices
        testCase.assertNumElements(reach_analyzer.verticesDirectlyReachableFromVertex(ray_ndx), 1);
        testCase.assertNumElements(reach_analyzer.conesDirectlyReachableFromVertex(ray_ndx), 1);

        next_vert = reach_analyzer.verticesDirectlyReachableFromVertex(ray_ndx);
        next_cone = reach_analyzer.conesDirectlyReachableFromVertex(ray_ndx);
        
        % Check vertex to vertex gains.
        [min_gain, max_gain] = reach_analyzer.gainsFromVertexToVertex(ray_ndx, next_vert);
        testCase.assertLessThanOrEqual(min_gain, 1.0);
        testCase.assertGreaterThanOrEqual(max_gain, 1.0);
        
        % Check vertex to cone gains.
        [min_gain, max_gain] = reach_analyzer.gainsFromVertexToCone(ray_ndx, next_cone);
        testCase.assertLessThanOrEqual(min_gain, 1.0);
        testCase.assertGreaterThanOrEqual(max_gain, 1.0);
      end

      % No vertices are reachable from origin.
      testCase.assertNumElements(reach_analyzer.verticesDirectlyReachableFromVertex(conical_partition.origin_index), 0);

      for cone_ndx = conical_partition.cone_indices
        testCase.assertNumElements(reach_analyzer.verticesDirectlyReachableFromCone(cone_ndx), 1);

        next_vert = reach_analyzer.verticesDirectlyReachableFromCone(cone_ndx);

        % Check cone to vertex gains.
        [min_gain, max_gain] = reach_analyzer.gainsFromConeToVertex(cone_ndx, next_vert);
        testCase.assertLessThanOrEqual(min_gain, 1.0);
        testCase.assertGreaterThanOrEqual(max_gain, 1.0);
      end


%       testCase.assertEqual(size(reach_analyzer.getEdgesFromVerticesToVertices(), 1), n_cones);
%       testCase.assertEqual(size(reach_analyzer.getEdgesFromConesToVertices(), 1), n_cones);
%       testCase.assertEqual(size(reach_analyzer.getEdgesFromVerticesToCones(), 1), n_cones);


      
%       testCase.assertNotEmpty(reach_analyzer.restricted_reachable_sets_from_unit_sphere_in_cone);
%       testCase.assertNotEmpty(reach_analyzer.can_flow_from_vertex_into_cone);
%       testCase.assertNotEmpty(reach_analyzer.can_flow_from_cone_to_vertex);
% 
%       cones_transition_graph = reach_analyzer.contractEdgesToBetweenCones(1:n_cones, 1:n_cones);
%       testCase.assertEqual(cones_transition_graph.numEdges(), n_cones^2);
% 
% 
%       % TODO: Test cycles
%       % cycles = reach_analyzer.getCycles();
%       % [cycle_min_gains, cycle_max_gains, cycles_nodes, cycles_edges] = reach_analyzer.getCycleGains();
%       % testCase.assertNumElements(cycles, 1);
%       % testCase.assertEqual(cycles{1}, conical_partition.ray_indices);
%       %
%       % % Since the flow is purely rotational, the gain of any solution is 1.0, but since we are approximating the system, the max and min gains we compute will be over-conservative, meaning:
%       % % `
%       % % `   cycle_min_gains ≤ 1.0 ≤ cycle_max_gains
%       % % `
%       % testCase.assertLessThanOrEqual(cycle_min_gains, 1.0);
%       % testCase.assertGreaterThanOrEqual(cycle_max_gains, 1.0);

      % testCase.verifyFail("Test case needs to be implemented.");
    end % End of function.
  end % End of test methods

end % End of test methods block
