function transformed_angles = mapAnglesByLinearMap(A, angles)
  arguments(Input)
    A (:, :) double {pwintz.validators.mustBeSquare};
    angles  (1, :) double;
  end % End of Input arguments block
  arguments(Output)
    transformed_angles (1, :) double;
  end % End of Output arguments block
  
  dim = size(A, 1);
  switch dim
    case 2
      vectors = A * pwintz.math.angle2UnitVector(angles);
      transformed_angles  = pwintz.math.atan2(vectors);
      transformed_angles = mod(transformed_angles, 2*pi);
    otherwise
      error("mapAnglesByLinearMap:UNEXPECTED_CASE", "Unexpected case: dim=%s.", dim);
  end
end