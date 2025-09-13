function state_cone_angles = generateFlowBasedConeAngles2D(options)
  arguments(Input)
    options.maxStateConeAngle      (1, 1) double = 2*pi / 10;
    options.maxDerivativeConeAngle (1, 1) double = 2*pi / 20;
    options.flowMapMatrix     (2, 2) double = -eye(2, 2); % "A_c" in \dot x = A_c x.
    options.flowSetStartAngle (1, 1) double =  0;
    options.flowSetEndAngle   (1, 1) double = pi;
    options.verbose logical = false;
  end % End of Input arguments block
  
  A_c = options.flowMapMatrix;
  
  % Create the baseline derivative vertex angles. We'll add more state angles to this list to account for the flow and jump sets.
  n_min_derivative_cones = ceil(2*pi / options.maxDerivativeConeAngle);
  derivative_cone_angles = pwintz.arrays.range("start", 0, "end", 2*pi, "n_values", n_min_derivative_cones, "includeEnd", false);

  state_cone_angles = mapAnglesByInvsLinearMap(A_c, derivative_cone_angles);
  
  % Add the angles from the flow set, jump set, and image of the jump.
  state_cone_angles = [state_cone_angles, options.flowSetStartAngle, options.flowSetEndAngle];
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