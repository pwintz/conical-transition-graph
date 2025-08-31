function v_normalized = nrv(v)
  % Compute the normalized radial vector (NRV) of v. 
  % If v is an array in the form v = [v1, v2, ..., vn], then the result is 
  %   v_normalized = [nrv(v1), nrv(v2), ..., nrv(vn)].
  arguments(Input)
    v (2, :) double;
  end
  arguments(Output)
    v_normalized (2, :) double;
  end
  
  nonzero_ndxs = vecnorm(v) > 1e-10;
  v_normalized = zeros(size(v));
  v_normalized(:, nonzero_ndxs) = v(:, nonzero_ndxs) ./ vecnorm(v(:, nonzero_ndxs));

  pwintz.assertions.assertAllAreMembers(...
    vecnorm(v_normalized), ...  % The norms of every column ...
    [0, 1], ...                 % ... must be either 0 or 1.
    tolerance=1e-6...
  );
end % end function