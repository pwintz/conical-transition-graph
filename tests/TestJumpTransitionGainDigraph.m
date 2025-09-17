classdef TestJumpTransitionGainDigraph < matlab.unittest.TestCase
% ` runtests TestJumpTransitionGainDigraph
  methods (Test)
    function test_one_mode(testCase)
      % ⋘────────── Setup ───────────⋙
      conical_partition = ConicalPartition.fromNumSlices(4);

      % quad_1_ndx = 1;
      % quad_2_ndx = 2;
      % quad_3_ndx = 3;
      % quad_4_ndx = 4;
      % left_half_ndxs = [2, 3]
      % right_half_ndxs = [1, 4]

      jump_set_cone_ndxs = 1:3;
      jump_map_matrix    = [1, 0; 0, -1]; % reflect across x-axis.
      
      % ⋘────────── Execute ─────────⋙
      jump_specs = JumpSpecification(jump_set_cone_ndxs, jump_map_matrix);
      jump_graph = JumpTransitionGainDigraph(conical_partition, jump_specs);
            
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(jump_graph.jump_set_cone_ndxs{1, 1}, 1:3);
      testCase.assertEqual(jump_graph.jump_set_image_cone_ndxs{1, 1}, [2, 3, 4]);% In left half plane. 
      mode_ndx = 1; 
      testCase.assertTrue(jump_graph.hasEdgeFromConeToCone(mode_ndx, 1, mode_ndx, 4));
      testCase.assertTrue(jump_graph.hasEdgeFromConeToCone(mode_ndx, 2, mode_ndx, 3));
      testCase.assertFalse(jump_graph.hasEdgeFromConeToCone(mode_ndx, 1, mode_ndx, 2));
      testCase.assertFalse(jump_graph.hasEdgeFromConeToCone(mode_ndx, 3, mode_ndx, 4));
      testCase.assertFalse(jump_graph.hasEdgeFromConeToCone(mode_ndx, 4, mode_ndx, 3));

      cycles = jump_graph.getCycles();
      testCase.assertNumElements(cycles, 1); % One cycle. 

      % testCase.assertEqual(cycles{1}, [2, 3])
    end % End of function.

    function test_two_modes(testCase)
      % ⋘────────── Setup ───────────⋙
      conical_partition_1 = ConicalPartition.fromNumSlices(4);
      conical_partition_2 = ConicalPartition.fromNumSlices(8);

      % quad_1_ndx = 1;
      % quad_2_ndx = 2;
      % quad_3_ndx = 3;
      % quad_4_ndx = 4;
      % left_half_ndxs = [2, 3]
      % right_half_ndxs = [1, 4]

      jump_set_cone_ndxs_1_to_1 = 1:3;
      jump_set_cone_ndxs_1_to_2 = 4;
      jump_set_cone_ndxs_2_to_1 = [4, 5];
      jump_set_cone_ndxs_2_to_2 = conical_partition_2.cone_indices;
      jump_map_matrix_1_to_1    = [1, 0; 0, -1];
      jump_map_matrix_1_to_2    = [1, 2; 2, -1];
      jump_map_matrix_2_to_1    = [1, 0; 1, -1];
      jump_map_matrix_2_to_2    = -eye(2);
      
      % ⋘────────── Execute ─────────⋙
      jump_specs_1_to_1 = JumpSpecification(jump_set_cone_ndxs_1_to_1, jump_map_matrix_1_to_1);
      jump_specs_1_to_2 = JumpSpecification(jump_set_cone_ndxs_1_to_2, jump_map_matrix_1_to_2);
      jump_specs_2_to_1 = JumpSpecification(jump_set_cone_ndxs_2_to_1, jump_map_matrix_2_to_1);
      jump_specs_2_to_2 = JumpSpecification(jump_set_cone_ndxs_2_to_2, jump_map_matrix_2_to_2);

      jump_specs = [
        jump_specs_1_to_1, jump_specs_1_to_2; 
        jump_specs_2_to_1, jump_specs_2_to_2;
      ];

      conical_partitions = [conical_partition_1, conical_partition_2];

      jump_graph = JumpTransitionGainDigraph(conical_partitions, jump_specs);
            
      % ⋘────────── Verify ──────────⋙
      testCase.assertEqual(jump_graph.jump_set_cone_ndxs{1, 1}, 1:3);
      testCase.assertEqual(jump_graph.jump_set_image_cone_ndxs{1, 1}, [2, 3, 4]);% In left half plane. 
      mode_ndx = 1; 
      testCase.assertTrue(jump_graph.hasEdgeFromConeToCone(mode_ndx, 1, mode_ndx, 4));
      testCase.assertTrue(jump_graph.hasEdgeFromConeToCone(mode_ndx, 2, mode_ndx, 3));
      testCase.assertFalse(jump_graph.hasEdgeFromConeToCone(mode_ndx, 1, mode_ndx, 2));
      testCase.assertFalse(jump_graph.hasEdgeFromConeToCone(mode_ndx, 3, mode_ndx, 4));
      testCase.assertFalse(jump_graph.hasEdgeFromConeToCone(mode_ndx, 4, mode_ndx, 3));


      mode_ndx = 2; 
      x = [1; 2];
      cone_ndx_1 = conical_partition_2.getConesContainingPoint( x);
      cone_ndx_2 = conical_partition_2.getConesContainingPoint(-x);

      testCase.assertTrue(jump_graph.hasEdgeFromConeToCone(mode_ndx, cone_ndx_1, mode_ndx, cone_ndx_2));
      testCase.assertTrue(jump_graph.hasEdgeFromConeToCone(mode_ndx, cone_ndx_2, mode_ndx, cone_ndx_1));
      testCase.assertFalse(jump_graph.hasEdgeFromConeToCone(mode_ndx, 1, mode_ndx, 1));
      testCase.assertFalse(jump_graph.hasEdgeFromConeToCone(mode_ndx, 1, mode_ndx, 2));

      % Check that getCycles doesn't produce errors.
      jump_graph.getCycles();
      % testCase.assertNumElements(cycles, 1); % One cycle. 

      % testCase.assertEqual(cycles{1}, [2, 3])
      
    end % End of function.
  end % End of test methods

end % End of test methods block
