classdef TestConicAbstraction < matlab.unittest.TestCase
  % ! You can run these tests using ConicAbstraction.test();
  % ! If you use the "runTestsStaticFunction" snippet to add a static method to ConicAbstraction.m

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

  % ╭───────────────────────────────────╮
  % │ ╭───────────────────────────────╮ │
  % │ │             Tests             │ │
  % │ ╰───────────────────────────────╯ │
  % ╰───────────────────────────────────╯
  methods (Test)
    function test_fromAngles_withoutWrappingPast2pi(testCase)
      % ⋘────────── Setup ───────────⋙
      A_c = eye(2);
      A_d = eye(2);
      flow_set_angles = [0, pi+0.1];
      jump_set_angles = [pi-0.1, pi+0.01];
      conic_abstraction = ConicAbstraction.fromAngles(n=12, flowMapMatrix = A_c, jumpMapMatrix = A_d, flowSetAngles=flow_set_angles, jumpSetAngles = jump_set_angles);
      % conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(4);
      % conic_abstraction = ConicAbstraction(conical_partition, );
      
      % ⋘────────── Check that the angles for the edges of the flow and jump sets have vertices─────────⋙
      testCase.assertTrue(ismembertol(flow_set_angles(1), conic_abstraction.conical_partition.state_vertex_angles));
      testCase.assertTrue(ismembertol(flow_set_angles(2), conic_abstraction.conical_partition.state_vertex_angles));
      testCase.assertTrue(ismembertol(jump_set_angles(1), conic_abstraction.conical_partition.state_vertex_angles));
      testCase.assertTrue(ismembertol(jump_set_angles(2), conic_abstraction.conical_partition.state_vertex_angles));
      
    end % End of function.
    
    function test_fromAngles_wrappingPast2pi(testCase)
      % ⋘────────── Setup ───────────⋙
      A_c =  eye(2);
      A_d = -eye(2);
      flow_set_angles = [-0.2, pi/2];
      jump_set_angles = [pi+0.1, 2*pi+0.01];

      % Compute the angles for the jump set image vertices.
      jump_set_image_vertices = A_d * pwintz.math.angle2UnitVector(jump_set_angles);
      jump_set_image_angles = pwintz.math.atan2(jump_set_image_vertices);

      conic_abstraction = ConicAbstraction.fromAngles(n=12, flowMapMatrix = A_c, jumpMapMatrix = A_d, flowSetAngles=flow_set_angles, jumpSetAngles = jump_set_angles);
      % conical_partition = ConicalPartition.fromNumberOfUniformDerivativeSlices(4);
      % conic_abstraction = ConicAbstraction(conical_partition, );
      
      % ⋘────────── Check that the angles for the edges of the flow and jump sets have vertices─────────⋙
      % isAngleInList = @(angle, angle_list) min(pwintz.math.angleDistance(angle, angle_list)) < ;
      
      state_vertex_angles = conic_abstraction.conical_partition.state_vertex_angles;
      test_cases = [ ... 
        struct("name", "flow_set_angles",       "angles", flow_set_angles), ...
        struct("name", "jump_set_angles",       "angles", jump_set_angles), ...
        struct("name", "jump_set_image_angles", "angles", jump_set_image_angles) ...
      ];
      for test_case = test_cases
        for i = 1:2
          angle = test_case.angles(i);
          angle_from_a_state_vertex_angle = min(pwintz.math.angleDistance(angle, state_vertex_angles));
          message = sprintf("%s(%d) mod 2pi = %f mod 2pi = %f is not in %s", test_case.name, i, angle, mod(angle, 2*pi), mat2str(state_vertex_angles));
          testCase.assertLessThan(angle_from_a_state_vertex_angle, ConicalPartition.ANGLE_TOL, message);
        end
      end

      % testCase.assertTrue(isAngleInList(flow_set_angles(1),       state_vertex_angles), ,       mat2str(state_vertex_angles)));
      % testCase.assertTrue(isAngleInList(flow_set_angles(2),       state_vertex_angles), sprintf("flow_set_angles(2) = %d is not in %s",      mod(flow_set_angles(2), 2*pi),       mat2str(state_vertex_angles)));
      % testCase.assertTrue(isAngleInList(jump_set_angles(1),       state_vertex_angles), sprintf("jump_set_angles(1) = %d is not in %s",      mod(jump_set_angles(1), 2*pi),       mat2str(state_vertex_angles)));
      % testCase.assertTrue(isAngleInList(jump_set_angles(2),       state_vertex_angles), sprintf("jump_set_angles(2) = %d is not in %s",      mod(jump_set_angles(2), 2*pi),       mat2str(state_vertex_angles)));
      % testCase.assertTrue(isAngleInList(jump_set_image_angles(1), state_vertex_angles), sprintf("jump_set_image_angles(1) = %d is not in %s",mod(jump_set_image_angles(1), 2*pi), mat2str(state_vertex_angles)));
      % testCase.assertTrue(isAngleInList(jump_set_image_angles(2), state_vertex_angles), sprintf("jump_set_image_angles(2) = %d is not in %s",mod(jump_set_image_angles(2), 2*pi), mat2str(state_vertex_angles)));
    end % End of function.
  end % End of test methods

end
