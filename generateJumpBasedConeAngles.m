function state_cone_angles = generateJumpBasedConeAngles(options)
  arguments(Input)
    options.maxStateConeAngle      (1, 1) double = 2*pi / 10;
    options.maxDerivativeConeAngle (1, 1) double = 2*pi / 20;
    options.jumpMapMatrix     (2, 2) double = -eye(2, 2);
    options.jumpSetStartAngle (1, 1) double =  0;
    options.jumpSetEndAngle   (1, 1) double = pi;
    options.verbose logical = false;
  end % End of Input arguments block
  
  A_d = options.jumpMapMatrix;
  
  % Compute the image of the jump set under the jump map.
  jump_set_boundary_angles = [options.jumpSetStartAngle, options.jumpSetEndAngle];

  % Create the baseline derivative vertex angles. We'll add more state angles to this list to account for the jump and jump sets.
  
  % Add the angles from the jump set, jump set, and image of the jump.
  state_cone_angles = [jump_set_boundary_angles];
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