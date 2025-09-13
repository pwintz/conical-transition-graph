function state_cone_angles = generateConeAngles2D(options)
  arguments(Input)
    options.maxStateConeAngle      (1, 1) double = 2*pi / 10;
    options.maxDerivativeConeAngle (1, 1) double = 2*pi / 20;
    options.flowMapMatrix     (2, 2) double = -eye(2, 2); % "A_c" in \dot x = A_c x.
    options.jumpMapMatrix     (2, 2) double =  eye(2, 2); % "A_d" in    x^+ = A_d x.
    options.flowSetStartAngle (1, 1) double =  0;
    options.flowSetEndAngle   (1, 1) double = pi;
    options.jumpSetStartAngle (1, 1) double = pi-0.1;
    options.jumpSetEndAngle   (1, 1) double = pi+0.1;
    options.verbose logical = false;
  end % End of Input arguments block
  
  A_c = options.flowMapMatrix;
  A_d = options.jumpMapMatrix;
  % ConicAbstraction.checkMatrices(A_c, A_d);
  % 
  % ConicAbstraction.checkSetAngles(options.flowSetAngles);
  % ConicAbstraction.checkSetAngles(options.jumpSetAngles);
  
  % Compute the image of the jump set under the jump map.
  jump_angles_range = [options.jumpSetStartAngle, options.jumpSetEndAngle];
  flow_angles_range = [options.flowSetStartAngle, options.flowSetEndAngle];

  jump_set_image_angles = mapAnglesByLinearMap(A_d, jump_angles_range);
  
  % Sort the two angles in jump_set_image_angles so that the angle from the first to the second is less than pi in the CCW direction. 
  if pwintz.math.angleDiffCCW(jump_set_image_angles, index=1) > pi
    jump_set_image_angles = circshift(jump_set_image_angles, 1);
  end
  assert(pwintz.math.angleDiffCCW(jump_set_image_angles, index=1) <= pi, ...
    'The jump_set_image_angles=%s should be ordered now such that the CCW distance from the first to second entry is < pi.')
  
  % Create the baseline derivative vertex angles. We'll add more state angles to this list to account for the flow and jump sets.
  n_min_derivative_cones = ceil(2*pi / options.maxDerivativeConeAngle);
  derivative_cone_angles = pwintz.arrays.range("start", 0, "end", 2*pi, "n_values", n_min_derivative_cones, "includeEnd", false);

  % derivative_cone_vertices = pwintz.math.angle2UnitVector(derivative_cone_angles);
  % state_cone_vertices = A_c \ derivative_cone_vertices;
  % state_cone_angles = pwintz.math.atan2(state_cone_vertices);
  state_cone_angles = mapAnglesByInvsLinearMap(A_c, derivative_cone_angles);
  
  % Add the angles from the flow set, jump set, and image of the jump.
  state_cone_angles = [state_cone_angles, flow_angles_range, jump_angles_range, jump_set_image_angles];
  % Update range to be [0, 2*pi)
  state_cone_angles = mod(state_cone_angles, 2*pi);
  % Sort and remove duplicates.
  state_cone_angles = unique(state_cone_angles, 'sorted');
  
  % ⋘──────── If any angles are too large, insert needed nodes ────────⋙
  angle_diff = pwintz.math.angleDiffCCW(state_cone_angles);
  extra_angles = [];
  for i = find(angle_diff > options.maxStateConeAngle)
    % We use a crude method for inserting more state cones. We simple take the last angle before the offending interval and step forward at the maximum allowed angle step. 
    extra_angles = [extra_angles, state_cone_angles(i) + (0:options.maxStateConeAngle:angle_diff(i))]; %#ok<AGROW>
  end
  state_cone_angles = unique([state_cone_angles, extra_angles], 'sorted');
  
end