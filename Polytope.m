
classdef Polytope < ConvexPolyhedron
  % Compact and convex polyhedron. Uses vertex representation or halfspace representation. 
 
   properties(SetAccess = immutable)
     % Define instance constants.
     vertices  (:, :) double;
     n_vertices; % Use this.count() to get number vertices
   end
 
  methods(Static)

    function polytope = getBox(dimension, options)
      arguments(Input)
        dimension (1,1) {mustBeInteger,mustBePositive};
        options.halfWidth (1, 1) double {mustBePositive} = 1.0;
      end % End of Input arguments block
      
      halfspace_representation = HalfspaceRepresentation.getBox(dimension, options.halfWidth);
      polytope = Polytope(halfspace_representation);
    end

    function convex_polyhedron = fromTriangle(points)
      dimension = size(points, 1);
      assert(ismember(dimension, [2, 3]), "Each column in 'points' must have 2 or 3 elements. Instead it had %d.", size(points, 1));
      assert(size(points, 2) == 3, "A triangle has three points. Instead %d were given.", size(points, 2))
      % TODO: Check colinear.
      
      convex_polyhedron = Polytope.fromConvexHull(points);
    end

    function convex_polyhedron = fromTetrahedron(points)
      dimension = size(points, 1);
      assert(dimension == 3, "Each column in 'points' must have three elements. Instead it had %d.", size(points, 1));
      assert(size(points, 2) == 4, "A tetrahedron has four points. Instead %d were given.", size(points, 2))
      % TODO: Check nonplanar.
      convex_polyhedron = Polytope.fromConvexHull(points);
    end

    function convex_polyhedron = fromConvexHull(points)
      dimension = size(points, 1);
      assert(ismember(dimension, [2, 3]), "Each column in 'points' must have 2 or 3 elements. Instead it had %d.", size(points, 1));
      % indices = convhull(points', 'simplify', true);
      % The output of convhull is an array of indices from the "points" list such that the points are arranged counter-clockwise, with the first one repeated at the end. We trim the last index so that 
      % vertices = points(:, indices(1:end-1));
      halfspace_representation = HalfspaceRepresentation.fromConvexHull(points);
      convex_polyhedron = Polytope(halfspace_representation);
    end
    
    % We define 
    function convex_polyhedron = fromPoint(point)
      assert(iscolumn(point), "point must be a column");
      assert(size(point, 1) == 2, "Each column in point must have two elements. Instead it had %d.", size(point, 1));
      halfspace_representation = HalfspaceRepresentation.fromConvexHull(point);
      convex_polyhedron = Polytope(halfspace_representation);
    end
    
    function convex_polyhedron = fromLineSegment(p1, p2)
      assert(iscolumn(p1), "p1 must be a column vector");
      assert(iscolumn(p2), "p2 must be a column vector");
      assert(norm(p1 - p2) > 0, "points must be distinct");
      convex_polyhedron = Polytope([p1, p2]);
    end
  end

   methods
    % Constructor
    function this = Polytope(halfspace_representation)
      arguments(Input)
        halfspace_representation (1,1) HalfspaceRepresentation {mustBeScalarOrEmpty} = [];
      end % End of Input arguments block
      

      if nargin() == 0 || isempty(halfspace_representation)
        halfspace_representation = HalfspaceRepresentation();
        vertices = double.empty(2, 0);
      else
        assert(halfspace_representation.isBounded(), "Half space representation of a polytope must define a bounded set.");
        vertices = halfspace_representation.getPolytopeVertices();
        assert(~isempty(vertices));
      end

      this = this@ConvexPolyhedron(halfspace_representation);

      % ⋘──────── Set the vertices properties ────────⋙
      
      this.vertices =  vertices;
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

    function [is_vertex, v_ndx] = isVertex(this, v)
      [is_vertex, v_ndx] = pwintz.arrays.isColumnIn(v, this.vertices, tolerance=1e-8);
      assert(numel(v_ndx) <= 1, 'Vertices are unique, so only one vertex expected.');
    end % End of function

     % ╭───────────────────────────────────────────────╮
     % │ ╭───────────────────────────────────────────╮ │
     % │ │             Vertex Operations             │ │
     % │ ╰───────────────────────────────────────────╯ │
     % ╰───────────────────────────────────────────────╯
     
     function poly = removeVertex(this, v)
       [is_member, v_ndx] = this.isVertex(v);
       if ~is_member
         warning("Polytope:notAVertex","The vertex %s was not in this polyhedron, which has vertices=%s.", mat2str(v), mat2str(this.vertices));
       end
       vertices_without_v  = this.vertices(:, setdiff(1:end, v_ndx));
       poly = ConvexPolyhedron(vertices_without_v);
     end % End of function
 
     % ╭────────────────────────────────────────────╮
     % │ ╭────────────────────────────────────────╮ │
     % │ │             Set Operations             │ │
     % │ ╰────────────────────────────────────────╯ │
     % ╰────────────────────────────────────────────╯
%      function result = intersectRay(this, v)
%        % v is a vector in R^2 (v ≠ 0) that defines a ray from the origin.
%  
%        assert(iscolumn(v), "v must be a column vector");
%        assert(size(v, 1) == 2);
%  
%        % ⋘──────────  Use the geo2d package to compute intersection. ──────────⋙
%        % Create a "ray" in the format used by geo2d, namely
%        % [x0, y0, dx, dy], where ray is [x0; y0] + t[dx; dy] for all t>0.
%        ray = [0 0 v'];
%  
%        % TODO: intersectRayPolygon will not throw an error if ray and the polygon are swapped. We should add some error checking to catch that. 
%        intersection_points = intersectRayPolygon(ray, this.vertices);
%        
%        result = ConvexPolyhedron(intersection_points);
%      end
 

 
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
         str = sprintf('Polytope with 0 vertices (empty).\n');
       elseif this.n_vertices == 1
         str = sprintf('Polytope with 1 vertex: \n%s\n', formattedDisplayText(this.vertices));
       else
         str = sprintf('Polytope with %d vertices: \n%s\n', this.n_vertices, formattedDisplayText(this.vertices));
       end
     end
 
     function disp(this)
       if ~isscalar(this)
         fprintf('%dx%d Polytope array\n', size(this, 1), size(this, 2))
         return
       end
       if this.n_vertices == 0
         fprintf("Polytope with 0 vertices (empty).\n");
       elseif this.n_vertices == 1
         fprintf("Polytope with 1 vertex: \n%s\n", formattedDisplayText(this.vertices));
       else
         fprintf("Polytope with %d vertices: \n%s\n", this.n_vertices, formattedDisplayText(this.vertices));
       end
     end
 
% 
%     function plot(this, varargin)
% 
%       vertices = this.halfspace_representation.getPolytopeVertices();
% 
%       switch this.dimension
%         case 2
%           patch('XData', this.vertices(1, :),'YData', this.vertices(2, :), varargin{:});
%         case 3
%           patch('XData', this.vertices(1, :),'YData', this.vertices(2, :), varargin{:});
%         otherwise
%           error("Unexpected case: %s.", this.dimension);
%       end
%     end % End of function


    function plot(this, varargin)
    
      DEFAULT_COLOR = 0.0 * [1, 1, 1];
      
      switch this.dimension
        case 2
          switch this.n_vertices
            case 1

              plot(this.vertices(1, :), this.vertices(2, :), "Marker", "x");
            case 2
              plot(this.vertices(1, :), this.vertices(2, :), "LineWidth", 5);
            otherwise
              patch("XData", this.vertices(1,:), "YData", this.vertices(2, :), "FaceColor", DEFAULT_COLOR, "FaceAlpha", 0.2, varargin{:});
          end
        case 3
          % K = convhull(x,y,z);
          x = this.vertices(1, :);
          y = this.vertices(2, :);
          z = this.vertices(3, :); 
          plot3(x, y, z);
          % trisurf(K,x,y,z,'Facecolor','cyan')

%           x
%           y
%           z
%           [k1,av1] = convhull(x,y,z, "Simplify", true);
%           k1
% 
%           trisurf(k1,x,y,z,'FaceColor','cyan')
        otherwise
          error("Unexpected case: %s.", this.dimension);
      end
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
 
 
 