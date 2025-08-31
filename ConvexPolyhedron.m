classdef ConvexPolyhedron < handle % & matlab.mixin.CustomDisplay
  % Tests for this class are located in TestConvexPolyhedron and can be run using ConvexPolyhedron.test().


  properties(SetAccess = immutable)
    % Define instance constants.
    vertices (2, :) double;
    dimension (1, 1) int32;
  end
  properties(SetAccess = immutable, GetAccess = private)
    % Define private variables.
    % first_vertex_ndx
    % center
    n_vertices; % Use this.count() to get number vertices
  end

  methods(Static)
    function convex_polyhedron = fromConvexHull(points)
      dimension = size(points, 1);
      assert(dimension == 2, "Each column in 'points' must have two elements. Instead it had %d.", size(points, 1));
      indices = convhull(points', 'simplify', true);
      % The output of convhull is an array of indices from the "points" list such that the points are arranged counter-clockwise, with the first one repeated at the end. We trim the last index so that 
      vertices = points(:, indices(1:end-1));
      convex_polyhedron = ConvexPolyhedron(vertices);
    end
    
    % We define 
    function convex_polyhedron = fromPoint(point)
      assert(iscolumn(point), "point must be a column");
      assert(size(point, 1) == 2, "Each column in point must have two elements. Instead it had %d.", size(point, 1));
      convex_polyhedron = ConvexPolyhedron(point);
    end
    
    function convex_polyhedron = fromLine(p1, p2)
      assert(iscolumn(p1), "p1 must be a column vector");
      assert(iscolumn(p2), "p2 must be a column vector");
      assert(norm(p1 - p2) > 0, "points must be distinct");
      convex_polyhedron = ConvexPolyhedron([p1, p2]);
    end

    function array = arrayOfEmpty(m, n)
      array = ConvexPolyhedron.empty(m, 0);
      array(1:m, 1:n) = ConvexPolyhedron();
    end

    % function convex_polyhedron = fromTriangle(p1, p2, p3)
    %   assert(iscolumn(p1), "p1 must be a column vector")
    %   assert(iscolumn(p2), "p2 must be a column vector")
    %   assert(iscolumn(p3), "p3 must be a column vector")
    %   assert(norm(p1 - p2) > 0, "points must be distinct")
    %   assert(norm(p1 - p2) > 0, "points must be distinct")
    %   convex_polyhedron = ConvexPolyhedron([p1, p2]);
    % end

    function tests(varargin) % Define convenience functions for running tests.
      TestConvexPolyhedron.runTests(varargin{:});
    end % End of function
  end

  methods(Access=private)
    % Constructor
    function this = ConvexPolyhedron(vertices)
      if nargin() == 0 || isempty(vertices)
        vertices = double.empty(2, 0);
      end

      dimension = size(vertices, 1);
      assert(dimension == 2, "Each column in vertices must have two elements. Instead size(vertices) = %s.", size(vertices));
      this.vertices =  vertices;
      this.dimension = dimension;
      % ! If isempty(vertices), then we can still have size(vertices, 2) > 0, since the shape of vertices can be 0xN for N > 0.
      if isempty(vertices)
        this.n_vertices = 0;
      else
        this.n_vertices = size(vertices, 2);
      end

      % TODO: Check that polyhedron is convex.
%       this.center = mean(vertices, 1);
% 
%       % ⋘────────── Find lexicographically "first" vertex ──────────⋙
%       [~, y_min_ndxs] = min(vertices(2,:));
%       if numel(y_min_ndxs) == 1
%         this.first_vertex_ndx = y_min_ndxs;
%       else
%         [~, x_min_ndx] = min(vertices(1, y_min_ndxs));
%         assert(numel(x_min_ndx) == 1, "There should be a unique 'first' vertex")
%         this.first_vertex_ndx = y_min_ndxs(x_min_ndx);
%       end
%       this.first_vertex_ndx 
    end
  end

  methods
    % ╭──────────────────────────────────────╮
    % │ ╭──────────────────────────────────╮ │
    % │ │             Queries              │ │
    % │ ╰──────────────────────────────────╯ │
    % ╰──────────────────────────────────────╯
    function v = getVertex(this, ndx)
      v = this.vertices(:, ndx);
    end % End of function

    function [is_member, v_ndx] = isVertex(this, v)
      [is_member, v_ndx] = ismembertol(v', this.vertices', 1e-6, 'ByRows', true);
      assert(size(v_ndx, 2) <= 1, 'Vertices are unique, show only one vertex expected.');
    end % End of function

    function TF = contains(this, points)
      % Use "isPointInPolygon" from the matGeom package.
      TF =  isPointInPolygon(points', this.vertices')';
    end % End of function

    % ╭───────────────────────────────────────────────╮
    % │ ╭───────────────────────────────────────────╮ │
    % │ │             Vertex Operations             │ │
    % │ ╰───────────────────────────────────────────╯ │
    % ╰───────────────────────────────────────────────╯
    
    function poly = Vertex(this, v)
      [is_member, v_ndx] = this.isVertex(v);
      if ~is_member
        warning("ConvexPolyhedron:notAVertex","The vertex %s was not in this polyhedron, which has vertices=%s.", mat2str(v), mat2str(this.vertices));
      end
      vertices_without_v  = this.vertices(:, setdiff(1:end, v_ndx));
      poly = ConvexPolyhedron(vertices_without_v);
    end % End of function

    % ╭────────────────────────────────────────────╮
    % │ ╭────────────────────────────────────────╮ │
    % │ │             Set Operations             │ │
    % │ ╰────────────────────────────────────────╯ │
    % ╰────────────────────────────────────────────╯
    function result = minkowskiSum(left, right)
      % ! We do a simple and inefficient implementation of the Minkowski sum here. 
      % ! It should be fine for small numbers of points, but the complexity is O(m*n) (or worse) instead of O(m + n), where m and n are the numbers of vertices in each polygon.

      % ⋘────────── Get Vertices of Operands ──────────⋙
      left_vertices   = ConvexPolyhedron.verticesFromOperand(left);
      right_vertices  = ConvexPolyhedron.verticesFromOperand(right);

      % ⋘────────── Get numbers of vertices ──────────⋙
      n_left = size(left_vertices, 2);
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
      % ⋘────────── Define a utility function ──────────⋙
      uniqueCols = @(array) unique(array', 'rows', 'stable')';
      
      % ⋘────────── Get vertices of operands ──────────⋙
      left_vertices   = ConvexPolyhedron.verticesFromOperand(left);
      right_vertices  = ConvexPolyhedron.verticesFromOperand(right);
      
      lastwarn(''); % Clear the last warning.
      if size(left_vertices, 2) == 2
        result = right.intersectRay(left_vertices);
      elseif size(right_vertices, 2) == 2
        result = left.intersectRay(right_vertices);
      else
        vertices1 = uniqueCols(left_vertices);
        vertices2 = uniqueCols(right_vertices);
        poly1 = polyshape(vertices1(1,:), vertices1(2,:));
        poly2 = polyshape(vertices2(1,:), vertices2(2,:));
        poly_intersection = intersect(poly1, poly2);
        intersection_vertices = poly_intersection.Vertices';
        result = ConvexPolyhedron(intersection_vertices);
      end
      warnmsg = lastwarn(); % Get the last warning
      if warnmsg ~= ""
        fprintf('There was a warning produced while finding the intersection between \n%sand \n%s\n', left, right);
      end
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
      assert(cond(A) < 1e9, "Only implemented for invertible A, for now. Otherwise we need to set collapses to a lower.");
      result = ConvexPolyhedron(A*this.vertices);
    end

    % ╭──────────────────────────────────────────────────────╮
    % │ ╭──────────────────────────────────────────────────╮ │
    % │ │             Overload Default Methods             │ │
    % │ ╰──────────────────────────────────────────────────╯ │
    % ╰──────────────────────────────────────────────────────╯
    % ╭───────────────────────────────────────────────────────────────────────╮
    % │             Make ConvexPolyhedron Behave Like a Container             │
    % ╰───────────────────────────────────────────────────────────────────────╯
    function result = isempty(this)
      result = isempty(this.vertices);
    end
    function result = count(this)
      result = size(this.vertices, 2);
    end

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
      left_subset_right = all(ismember( left.vertices', right.vertices', 'rows'));
      right_subset_left = all(ismember(right.vertices',  left.vertices', 'rows'));
      is_equal = left_subset_right && right_subset_left;
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
      
      if this.n_vertices == 0
        str = sprintf('ConvexPolyhedron with 0 vertices (empty).\n');
      elseif this.n_vertices == 1
        str = sprintf('ConvexPolyhedron with 1 vertex: \n%s\n', formattedDisplayText(this.vertices));
      else
        str = sprintf('ConvexPolyhedron with %d vertices: \n%s\n', this.n_vertices, formattedDisplayText(this.vertices));
      end
    end

    function disp(this)
      if ~isscalar(this)
        fprintf('%dx%d ConvexPolyhedron array\n', size(this, 1), size(this, 2))
        return
      end
      if this.n_vertices == 0
        fprintf("ConvexPolyhedron with 0 vertices (empty).\n");
      elseif this.n_vertices == 1
        fprintf("ConvexPolyhedron with 1 vertex: \n%s\n", formattedDisplayText(this.vertices));
      else
        fprintf("ConvexPolyhedron with %d vertices: \n%s\n", this.n_vertices, formattedDisplayText(this.vertices));
      end
    end

    function out = isequaln(left, right)
      out = eq(left, right);
    end % End of function

    function plot(this, varargin)
      patch('XData', this.vertices(1, :),'YData', this.vertices(2, :), varargin{:});
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


