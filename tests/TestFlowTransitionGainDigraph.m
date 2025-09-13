classdef TestFlowTransitionGainDigraph < matlab.unittest.TestCase
% ` runtests TestFlowTransitionGainDigraph

  methods (Test)
    function test_constructor_2D(testCase)
      % ⋘────────── Setup ───────────⋙
      conical_partition = ConicalPartition.fromNumSlices(4);
      flow_set_cone_ndxs = 1:(conical_partition.n_cones / 2)
      flow_map_matrix    = magic(2);   
      
      % ⋘────────── Execute ─────────⋙
      flow_digraph = FlowTransitionGainDigraph(conical_partition, flow_set_cone_ndxs, flow_map_matrix);
      
      
      % ⋘────────── Verify ──────────⋙
      % testCase.verifyFail("Test case needs to be implemented.");
    end % End of function.
  end % End of test methods

end % End of test methods block
