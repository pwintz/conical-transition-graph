classdef ConvexPolyhedron < handle
  % Tests for this class are located in TestConvexPolyhedron and can be run using 
  % ` runtests TestConvexPolyhedron

  properties(SetAccess = immutable, GetAccess = protected)
    % Define instance constants.
    dimension (1, 1) int32;
    halfspace_representation (1, 1) HalfspaceRepresentation;
  end

  methods(Static)

    function array = arrayOfEmpty(m, n)
      array = ConvexPolyhedron.empty(m, 0);
      array(1:m, 1:n) = ConvexPolyhedron();
    end

  end

  methods
    % Constructor
    function this = ConvexPolyhedron(halfspace_representation)
      if nargin() == 0
        halfspace_representation = HalfspaceRepresentation();
      end
      dimension = halfspace_representation.ambient_dimension;
      % assert(ismember(dimension, [2,3]), "Each column in vertices must have two or three elements. Instead dimension = %d.", dimension);
      this.halfspace_representation = halfspace_representation;
      this.dimension = dimension;
    end
  end

  methods
    % ╭──────────────────────────────────────╮
    % │ ╭──────────────────────────────────╮ │
    % │ │             Queries              │ │
    % │ ╰──────────────────────────────────╯ │
    % ╰──────────────────────────────────────╯

    function is_point_in_poly = contains(this, points)
      % Use "isPointInPolygon" from the matGeom package.
      is_point_in_poly = this.halfspace_representation.containsPoints(points);
      % TF =  isPointInPolygon(points', this.vertices')';
    end % End of function

    function is_bounded = isBounded(this)
      is_bounded = this.halfspace_representation.isBounded();
    end

    function is_unbounded = isUnbounded(this)
      is_unbounded = ~this.halfspace_representation.isBounded();
    end

    % ╭────────────────────────────────────────────╮
    % │ ╭────────────────────────────────────────╮ │
    % │ │             Set Operations             │ │
    % │ ╰────────────────────────────────────────╯ │
    % ╰────────────────────────────────────────────╯
    function result = minkowskiSum(left, right)
      % ! We do a simple and inefficient implementation of the Minkowski sum here. 
      % ! It should be fine for small numbers of points, but the complexity is O(m*n) (or worse) instead of O(m + n), where m and n are the numbers of vertices in each polygon.

      if isa(right, "ConvexPolyhedralCone")
        result = right.minkowskiSum(left);
        return
      end
      disp("Using Polytope/minkowskiSum for sum of two polytopes.");
      error("Polytope/minkowskiSum is not implemented.");

      % ⋘────────── Get Vertices of Operands ──────────⋙
      left_vertices   = ConvexPolyhedron.verticesFromOperand(left);
      right_vertices  = ConvexPolyhedron.verticesFromOperand(right);

      % ⋘────────── Get numbers of vertices ──────────⋙
      n_left  = size(left_vertices, 2);
      n_right = size(right_vertices, 2);
      n_points = n_left * n_right;

      % ⋘────────── Generate all the combinations of sums ──────────⋙
      points = nan(2, n_points);
      i_point = 1;
      for i_left = 1:n_left
        for i_right = 1:n_right
          points(:, i_point) = left_vertices(:, i_left) + right_vertices(:, i_right);
          i_point = i_point + 1;
        end
      end

      % ⋘────────── Set result to convex hull ──────────⋙
      result = ConvexPolyhedron.fromConvexHull(points);

      % # The following code was an attempt at a more efficient implementation.
      % https://cp-algorithms.com/geometry/minkowski.html
      % left_ndx  = left.first_vertex_ndx;
      % right_ndx = right.first_vertex_ndx;
      % n_result_vertices = left.n_vertices + right.n_vertices;
      % for result_ndx = 1:n_result_vertices
      %   result_vertices(result_ndx) = left.vertices(left_ndx) + right.vertices(right_ndx);
      %   if atan2(y, x)
      % end
    end

    function result = intersection(left, right)
      intersection_halfspace_representation = left.halfspace_representation.intersection(right.halfspace_representation);
      if intersection_halfspace_representation.isBounded()
        result = Polytope(intersection_halfspace_representation);
      else
        result = ConvexPolyhedron(intersection_halfspace_representation);
      end

      return
      

      % % ⋘────────── Get vertices of operands ──────────⋙
      % left_vertices   = ConvexPolyhedron.verticesFromOperand(left);
      % right_vertices  = ConvexPolyhedron.verticesFromOperand(right);
      % 
      % lastwarn(''); % Clear the last warning.
      % if size(left_vertices, 2) == 2
      %   result = right.intersectRay(left_vertices);
      % elseif size(right_vertices, 2) == 2
      %   result = left.intersectRay(right_vertices);
      % else
      %   vertices1 = pwintz.arrays.uniqueColumns(left_vertices);
      %   vertices2 = pwintz.arrays.uniqueColumns(right_vertices);
      %   poly1 = polyshape(vertices1(1,:), vertices1(2,:));
      %   poly2 = polyshape(vertices2(1,:), vertices2(2,:));
      %   poly_intersection = intersect(poly1, poly2);
      %   intersection_vertices = poly_intersection.Vertices';
      %   result = ConvexPolyhedron(intersection_vertices);
      % end
      % warnmsg = lastwarn(); % Get the last warning
      % if warnmsg ~= ""
      %   fprintf('There was a warning produced while finding the intersection between \n%sand \n%s\n', left, right);
      % end
    end

    function result = intersectRay(this, v)
      % v is a vector in R^2 (v ≠ 0) that defines a ray from the origin.

      assert(iscolumn(v), "v must be a column vector");
      assert(size(v, 1) == 2);

      % ⋘──────────  Use the geo2d package to compute intersection. ──────────⋙
      % Create a "ray" in the format used by geo2d, namely
      % [x0, y0, dx, dy], where ray is [x0; y0] + t[dx; dy] for all t>0.
      ray = [0 0 v'];

      % TODO: intersectRayPolygon will not throw an error if ray and the polygon are swapped. We should add some error checking to catch that. 
      intersection_points = intersectRayPolygon(ray, this.vertices);
      
      result = ConvexPolyhedron(intersection_points);
    end

    function result = linearTransform(this, A)
      error("Deprectated. Use a polytope or convexpolyhedralcone.");
      % assert(cond(A) < 1e9, "Only implemented for invertible A, for now. Otherwise we need to set collapses to a lower.");
      % result = ConvexPolyhedron(A*this.vertices);
    end

    function result = atXDoesVPointInward(this, bnd_point, vector)
      arguments(Input)
        this;
        bnd_point (:, 1) double;
        vector    (:, 1) double;
      end % End of Input arguments block
      arguments(Output)
        result logical {pwintz.validators.mustBeScalar};
      end % End of Output arguments block
      
      [A_active, ~] = this.halfspace_representation.activeInequalityConstraintsAtPoint(bnd_point);
      A_eq = this.halfspace_representation.A_eq;

      % The point "bnd_point" is on the boundary, so A_active*bnd_point = b_active.
      % We want to test if A_active*(bnd_point + vector) is <= or > b_active, which is the same as asking if 
      % A_active*vector) is positive or negative.
      result = all(A_active * vector <= 0) && all(abs(A_eq * vector) < 1e-7);

      % assert(this.halfspace_representation.n_equality_constraints == 0, "Not implemented for equality constraints.");
      % A_eq = this.halfspace_representation.A_eq;
      % assert(all(abs(A_eq * v) < 1e-6));
    end

    % ╭──────────────────────────────────────────────────────╮
    % │ ╭──────────────────────────────────────────────────╮ │
    % │ │             Overload Default Methods             │ │
    % │ ╰──────────────────────────────────────────────────╯ │
    % ╰──────────────────────────────────────────────────────╯
    % % ╭───────────────────────────────────────────────────────────────────────╮
    % % │             Make ConvexPolyhedron Behave Like a Container             │
    % % ╰───────────────────────────────────────────────────────────────────────╯
    % function result = isempty(this)
    %   result = isempty(this.vertices);
    % end
    % function result = count(this)
    %   result = size(this.vertices, 2);
    % end

    % ╭────────────────────────────────────────────╮
    % │             Overload Operators             │
    % ╰────────────────────────────────────────────╯
    % Overload "|" operator to compute the intersection of sets.
    function result = or(left, right)
      result = intersection(left, right);
    end

    % Overload "+" operator to compute the Minkowski sum of sets.
    function result = plus(left, right)
      result = minkowskiSum(left, right);
    end

    % Overload "*" operator
    function result = mtimes(A, convex_polyhedron)
      result = convex_polyhedron.linearTransform(A);
    end

    % Overload "==" operator
    function is_equal = eq(left, right)

      % These checks to not test allow for tolerance and depend on vertices, which is only available for polytopes.
      left_subset_right = all(ismember( left.vertices', right.vertices', 'rows'));
      right_subset_left = all(ismember(right.vertices',  left.vertices', 'rows'));
      is_equal = left_subset_right && right_subset_left;


      % We should use an equality check in the halfspace representation. 
      is_equal = left.halfspace_representation == right.halfspace_representation;


    end % End of function
    
    % Overload "~=" operator
    function is_not_equal = ne(left, right)
      is_not_equal = ~(left == right);
    end % End of function

    % ╭─────────────────────────────────────────────────────╮
    % │             Overload Built-in Functions             │
    % ╰─────────────────────────────────────────────────────╯
    function str = char(this)
      % This method defines the string representation that is 
      % inserted when using "%s" in a string format, such as 
      % fprintf("%s", polyhedron) or sprintf("%s", polyhedron); 
      % It does not change the "disp" output. 
      arguments(Output)
        str char; % Ensure output is cast to char, even if you create a string.
      end 
      
      str = sprintf('ConvexPolyhedron with halfspace representation: \n%s\n', this.halfspace_representation);
    end

    function disp(this)
      if ~isscalar(this)
        fprintf('%dx%d ConvexPolyhedron array\n', size(this, 1), size(this, 2))
        return
      end

      fprintf("%s\n", char(this))
    end

    function out = isequaln(left, right)
      out = eq(left, right);
    end % End of function
    
  end% End of methods block

  methods(Access=private, Static)
    function vertices = verticesFromOperand(operand)
      % Given an operand that is either a ConvexPolyhedron or a 2xN array, return the corresponding set of vertices. For the case where "operand" is an array, the "vertices" are the "operand" array.
      
      if isa(operand, "ConvexPolyhedron")
        vertices = operand.vertices;
      elseif isnumeric(operand) && size(operand, 1) == 2
        vertices = operand;
      else
        error("ConvexPolyhedron:invalidOperand","operand is invalid type (%s) or size (%s). Expected a Convex Polyhedron or a 2xN array.", class(operand), mat2str(size(operand)));
      end
    end
  end % End of private static methods block

  % methods(Access=protected)
  %   function groups = getPropertyGroups(obj)
  %       if isscalar(obj)
  %         % Specify the values to be displayed for properties
  %         propList = struct('vertices', newline + formattedDisplayText(obj.vertices));%) mat2str(obj.vertices, 3));
  %         groups = matlab.mixin.util.PropertyGroup(propList);
  %       else
  %         % Nonscalar case: call superclass method
  %         groups = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
  %       end
  %   end
  % end

  % TODO: Try to make nicer formatting for the strings.
  methods%(Access=private)
    % function str = verticesToString(this)
    %   string_array = strings(2, this.n_vertices)
    % end
  end
end


