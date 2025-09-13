function points = sphereGrid3d(options)
  % Create a grid of points that covers the unit sphere.
  arguments(Input)
    options.nLinesOfLongitude = 100;
    options.nLinesOfLatitude  = 100;
  end % End of Input arguments block
  

  % ⋘──────── Generate vertices from angles ────────⋙
  phis   = pwintz.arrays.range(0, 2*pi, n_values=options.nLinesOfLongitude, includeEnd=false);
  thetas = pwintz.arrays.range(0,   pi, n_values=options.nLinesOfLatitude,  includeEnd=true);

  [Phi, Theta] = meshgrid(phis, thetas);
  Phi   = reshape(Phi,   1, []); % Reshape into column vectors.
  Theta = reshape(Theta, 1, []); % Reshape into column vectors.

  spherical2cartesian = @(r, theta, phi) [
    r .* sin(theta) .* cos(phi); % x
    r .* sin(theta) .* sin(phi); % y
    r .* cos(theta);             % z
  ];

  points = spherical2cartesian(1, Theta, Phi);
  % points = pwintz.arrays.mapRows(@(theta, phi) spherical2cartesian(1, theta, phi)', Theta, Phi);

  pwintz.assertions.assertNumRows(points, 3);
  pwintz.assertions.assertNumColumns(points, options.nLinesOfLongitude * options.nLinesOfLatitude);

  % Check that all of the rows have norm 1. 
  pwintz.assertions.assertAlmostEqual(vecnorm(points), 1);

end % end function

