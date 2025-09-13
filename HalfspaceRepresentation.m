classdef HalfspaceRepresentation

  % ` runtests TestHalfspaceRepresentation
  
  properties % Define instance variables.
    A_ineq (:, :) double {mustHaveOnlyNonzeroRows};
    b_ineq (:, 1) double;
    A_eq   (:, :) double {mustHaveOnlyNonzeroRows};
    b_eq   (:, 1) double;

    % Derived quantities
    ambient_dimension        int32;
    n_inequality_constraints int32;
    n_equality_constraints   int32;

    % For unbounded sets, we use the following bounding box to numerically bound sets.
    bounding_box_half_width = 1e1;
  end % End of properties block

  properties(SetAccess = immutable, GetAccess = private) % Define private variables.
    vertices; % If a polytope, store the vertices.
    rays;     % If a polyhedral cone, store the rays.
    
    A_bounding_box;
    b_bounding_box;
    bounding_box_indices;
  end % End of private properties block

  properties(Dependent)
    % Define dependent variables.
    is_unconstrained;
  end
  methods
    function value = get.is_unconstrained(this)
      value = isempty(this.A_ineq) && isempty(this.A_eq);
    end
  end % End of dependent properties block

  methods(Static)

    function halfspace_rep = fromConvexHull(vertices)
      arguments(Input)
        vertices (:, :) double {mustBeFinite,mustBeNonNan};
      end % End of Input arguments block
      
      dimension = size(vertices, 1);
      n_vertices = size(vertices, 2);

      assert(ismember(dimension, [2,3]), pwintz.strings.format("Expected vertices to be dimension 2 or 3. Instead they had dimension %d. The size of vertices is %z", dimension, vertices));
      assert(n_vertices > 0, "half space representation should have at least one vertex. vertices=%s", mat2str(vertices))

      % Add the origin and then remove any duplicate vertices. 
      % vertices = pwintz.arvertices.uniqueColumns([vertices, origin], tolerance=1e-8);
    
      % Find the matrices and vectors defining the polyhedron.
      [A_ineq, b_ineq, A_eq, b_eq] = polyhedron.vert2lcon(vertices', 1e-6);
    
      
      halfspace_rep = HalfspaceRepresentation(A_ineq, b_ineq, A_eq, b_eq, vertices=vertices);
      % figure(1);
      % clf();
      % xlim(3*[-1, 1]);
      % ylim(3*[-1, 1]);
      % axis square;
      % hold on;
      % halfspace_rep.plot();
      assert(halfspace_rep.isBounded(), "A convex hull must be bounded.");
    end % End of function

    function halfspace_rep = fromConicalHull(rays)

      % Add the origin and then remove any duplicate rays. 
      origin = 0*rays(:, 1);
      rays = pwintz.arrays.uniqueColumns([rays, origin], tolerance=1e-8);

      % Find the matrices and vectors defining the polyhedron.
      [A_ineq, b_ineq, A_eq, b_eq] = polyhedron.vert2lcon(rays', 1e-6);

      is_conal_ineq_constraint = abs(b_ineq) < 1e-6;
      % ! Previously, we checked if there were multiple inequality constraints not going through the origin under the mistaken belief that there couldn't be more than one.  In fact, the bounding box can introduce multiple facets if the cone intersects an edge or corner of the box.

      A_ineq = A_ineq(is_conal_ineq_constraint, :); 
      b_ineq = b_ineq(is_conal_ineq_constraint);
      halfspace_rep = HalfspaceRepresentation(A_ineq, b_ineq, A_eq, b_eq, rays=rays);
    end % End of function

    function halfspace_representation = getSquare(half_width)
      arguments(Input)
        half_width = 1;
      end % End of Input arguments block
      
      dimension  = 2;
      A_ineq = [eye(dimension); -eye(dimension)];
      b_ineq = half_width * ones(2*dimension, 1);
      halfspace_representation = HalfspaceRepresentation(A_ineq, b_ineq);
    end
    

    function halfspace_representation = getBox(dimension, half_width)

      arguments(Input)
        dimension  = 3;
        half_width = 1;
      end % End of Input arguments block
      
      A_ineq = [eye(dimension); -eye(dimension)];
      b_ineq = half_width * ones(2*dimension, 1);
      halfspace_representation = HalfspaceRepresentation(A_ineq, b_ineq);
    end

  end % End static methods block

  methods
    
    % Constructor
    function this = HalfspaceRepresentation(A_ineq, b_ineq, A_eq, b_eq, options)
      arguments(Input)
        A_ineq (:, :) double {mustHaveOnlyNonzeroRows} = double.empty(0, 1);
        b_ineq (:, 1) double                           = double.empty(0, 1);
        A_eq   (:, :) double {mustHaveOnlyNonzeroRows} = double.empty(0, 1);
        b_eq   (:, 1) double                           = double.empty(0, 1);

        options.vertices = []; % For storing the vertices of a polytope.
        options.rays     = []; % For storing the rays defining a polyhedral cone. 
      end % End of Input arguments block
      
      % Check the number of arguments.
      if ~ismember(nargin(), [0, 2, 4])
        error("HalfspaceRepresentation. Unexpected number of arguments. Must be 0, 2, or 4. Instead was %d", nargin());
      end
      
      % Check that A_ineq and b_ineq are compatible
      pwintz.assertions.assertSameNumRows(this.A_ineq, this.b_ineq, leftName="A_ineq", rightName="b_ineq");
      % Check that A_eq and b_eq are compatible
      pwintz.assertions.assertSameNumRows(this.A_eq,   this.b_eq, leftName="A_eq", rightName="b_eq");

      if isempty(A_ineq)
        % If there are no inequality constraints, then make A_ineq have the right width, while having zero rows.
        A_ineq = double.empty(0, size(A_eq, 2));
      end

      ambient_dimension        = size(A_ineq, 2);

      % If the A_eq and b_eq arguments are not given
      if isempty(A_eq)
        A_eq = double.empty(0, ambient_dimension);
      end
      
      % ⋘──────── Check argument sizes ────────⋙
      [A_ineq, b_ineq, A_eq, b_eq] = preprocessLinearConstraints(A_ineq, b_ineq, A_eq, b_eq);

      n_inequality_constraints = size(A_ineq, 1);
      n_equality_constraints   = size(A_eq,   1);
      pwintz.assertions.assertSize(A_eq,   [  n_equality_constraints, ambient_dimension]);
      pwintz.assertions.assertSize(b_ineq, [n_inequality_constraints, 1]);
      pwintz.assertions.assertSize(b_eq,   [  n_equality_constraints, 1]);

      this.A_ineq = A_ineq;
      this.b_ineq = b_ineq;
      this.A_eq = A_eq;
      this.b_eq = b_eq;

% ! I don't thing the following lines are needed anymore.
%       Ab_ineq = [A_ineq, b_ineq];
%       Ab_eq   = [A_eq, b_eq];
% 
%       % We delete A_ineq b_ineq A_eq b_eq to prevent us from accidentally using them instead of Ab_ineq and Ab_eq.
%       clear A_ineq b_ineq A_eq b_eq; 
% 
%       % Remove duplicate rows
%       Ab_ineq = pwintz.arrays.uniqueRows(Ab_ineq, tolerance=1e-9);
% 
%       % Find any rows that where the negation of the row is present in the array.
%       [~, duplicate_row_ndxs, ~, unique_row_ndxs] = pwintz.arrays.duplicatedRows([Ab_ineq; -Ab_ineq], tolerance=1e-9);
% 
% 
%       % all_indices = sort([duplicate_row_ndxs; unique_row_ndxs])
% 
%       % unique_row_ndxs    = intersect(1:size(A_ineq, 1), unique_row_ndxs)
%       % duplicate_row_ndxs = intersect(1:size(A_ineq, 1), duplicate_row_ndxs)
% 
%       ineq_rows_to_add_to_eq_constraints    = intersect(1:size(Ab_ineq, 1), duplicate_row_ndxs);
%       ineq_rows_to_keep_in_ineq_constraints = intersect(1:size(Ab_ineq, 1), unique_row_ndxs);
%       
%       % Add rows from Ab_ineq that create equality constraints.
%       Ab_eq = [Ab_eq; Ab_ineq(ineq_rows_to_add_to_eq_constraints, :)];
% 
      % % ╭────────────────────────────────────────────────────────────────────╮
      % % │             Remove duplicate equality constraint rows.             │
      % % ╰────────────────────────────────────────────────────────────────────╯
      % % TODO: This will not remove redundant rows if they have opposite signs. 
      % Ab_eq = pwintz.arrays.uniqueRows(Ab_eq, tolerance=1e-9);
      % 
      % % Remove any rows from Ab_ineq that we didnt add 
      % Ab_ineq = Ab_ineq(ineq_rows_to_keep_in_ineq_constraints, :);

      % [Ab_ineq_unique, src_ndxs_cell, out_rows_ndxs_cell] = pwintz.arrays.uniqueRows(C)
% 
%       this.A_ineq = Ab_ineq(:, 1:end-1);
%       this.b_ineq = Ab_ineq(:, end);
%       this.A_eq   = Ab_eq(:, 1:end-1);
%       this.b_eq   = Ab_eq(:, end);


      % Recompute number of constraints.
      n_inequality_constraints = size(this.A_ineq, 1);
      n_equality_constraints   = size(this.A_eq,   1);

      this.ambient_dimension        = ambient_dimension;       
      this.n_inequality_constraints = n_inequality_constraints;
      this.n_equality_constraints   = n_equality_constraints;  


      % Define box constraints -b <= x <= b, as 
      % *   I*x <= b  
      % *  -I*x <= b  
      A_bounding_box = [eye(this.ambient_dimension); -eye(this.ambient_dimension)];
      b_bounding_box = this.bounding_box_half_width * ones(2*this.ambient_dimension, 1);

      % The indices of the bounding box constraints when appended after the inequality constraints.
      this.bounding_box_indices = this.n_inequality_constraints + (1:2*this.ambient_dimension);
      
      % Check that bounding box constraints are OK (include region around the origin).
      x_inside = ones(this.ambient_dimension, 1);
      assert(all(A_bounding_box *  x_inside <= b_bounding_box));
      assert(all(A_bounding_box * -x_inside <= b_bounding_box));
      x_outside = 2*this.bounding_box_half_width*x_inside;
      assert(any(A_bounding_box *  x_outside > b_bounding_box));
      assert(any(A_bounding_box * -x_outside > b_bounding_box));

      % Assign property values.
      this.A_bounding_box = A_bounding_box;
      this.b_bounding_box = b_bounding_box;

      this.vertices = options.vertices;
      this.rays     = options.rays;

      % Sanity check.
      pwintz.assertions.assertNumColumns(this.A_ineq, ambient_dimension, name="A_ineq");
      pwintz.assertions.assertNumColumns(this.A_eq,   ambient_dimension, name="A_eq");
      pwintz.assertions.assertNumRows(this.A_ineq, n_inequality_constraints, name="A_ineq");
      pwintz.assertions.assertNumRows(this.A_eq,   n_equality_constraints,   name="A_eq");
    end % End of constructor

    function vertices = getPolytopeVertices(this)

      [vertices, is_vertex_on_boundary_of_box] = getBoxBoundedVertices(this);
      if any(is_vertex_on_boundary_of_box)
        error("getPolytopeVertices(): This set is unbounded (or, at least, extends to the bounding box), so it is not a polytope.");
      end

%       A_with_bnd_box = [this.A_ineq; this.A_bounding_box];
%       b_with_bnd_box = [this.b_ineq; this.b_bounding_box];
% 
%       [vertices, not_redun_ineq] = polyhedron.lcon2vert(A_with_bnd_box, b_with_bnd_box, this.A_eq, this.b_eq);
%       
%       is_a_bounding_box_constraint_active = any(ismember(not_redun_ineq, this.bounding_box_indices));
%       if is_a_bounding_box_constraint_active
%         error("getPolytopeVertices(): This set is unbounded, so it is not a polytope.");
%       end

    end

    function rays = getConeRays(this)
      % If this is a cone, then find all of the rays from the origin. Otherwise, throw an error.
      % The rays are returned as an array with a ray in each column. The rays are normalized to unit length.

      if this.is_unconstrained
        rays = [];
        return
      end

      % Find the vertices of the polyhedron intersected with the (large) bounding box.
      [verts, is_vertex_on_boundary_of_box] = this.getBoxBoundedVertices();
      
      rays   = verts(:, is_vertex_on_boundary_of_box);
      origin = verts(:, ~is_vertex_on_boundary_of_box);

      pwintz.assertions.assertIsColumn(origin, "This set must be a cone, meaning there is only one bounded vertex, which is the origin.");

      % ⋘──────── Normalize the rays ────────⋙
      % Divide each column by the norm of that column.
      rays = rays ./ vecnorm(rays);
    end

    function is_in = containsPoint(this, point)
      assert(iscolumn(point), 'Expected a single column vector.');
      pwintz.assertions.assertNumRows(point, this.ambient_dimension, "point does not match the dimension of this halfspace representation.");
      tolerance = 1e-9;
      satisfies_ineq = all(this.A_ineq * point <= this.b_ineq + tolerance);
      satisfies_eq   = all(abs(this.A_eq * point - this.b_eq) < tolerance);

      is_in = satisfies_ineq && satisfies_eq;
    end

    function is_in = containsPoints(this, points)
      pwintz.assertions.assertNumRows(points, this.ambient_dimension, "The size of points do not match the dimension of this halfspace representation.");
      tolerance = 1e-9;
      satisfies_ineq = pwintz.logical.allPerColumn(this.A_ineq * points <= this.b_ineq + tolerance);
      satisfies_eq   = pwintz.logical.allPerColumn(abs(this.A_eq * points - this.b_eq) < tolerance);
      
      is_in = satisfies_ineq & satisfies_eq;
    end

    function [A, b] = activeInequalityConstraintsAtPoint(this, x)
      % Get the set of inequality constraints from "A_ineq * x < b_ineq" such that "A_ineq * x = b_ineq"  (up to numerical tolerance) at the given point x.
      TOL = 1e-7;
      assert(all(this.A_ineq * x <= this.b_ineq + TOL), "The point x=%s is not in this halfspace:\n%s\nWe have this.A_ineq * x = %s not less than this.b_ineq = %s.", mat2str(x), this, mat2str(this.A_ineq * x), mat2str(this.b_ineq));

      slack = this.b_ineq - this.A_ineq * x; % >= 0.
      assert(all(slack >= -TOL)); % All entries should be nonnegative.
      
      active_indices = slack < TOL;
      A = this.A_ineq(active_indices, :);
      b = this.b_ineq(active_indices);
    end % End of function

    function is_bounded = isBounded(this)
      if this.n_equality_constraints == 0 && this.n_inequality_constraints <= 1
        is_bounded = false;
        return;
      end

      % Find the vertices of the polyhedron intersected with the (large) bounding box.
      [~, is_vertex_on_boundary_of_box] = this.getBoxBoundedVertices();

      is_bounded = all(~is_vertex_on_boundary_of_box);
    end % End of function

    function result = intersection(this, other)
      % this_Ab = [this.A_ineq, this.b_ineq];
      % other_Ab = [other.A_ineq, other.b_ineq];
      A_ineq_both = [this.A_ineq; other.A_ineq];
      b_ineq_both = [this.b_ineq; other.b_ineq];
      A_eq_both   = [this.A_eq;   other.A_eq];
      b_eq_both   = [this.b_eq;   other.b_eq];

      % [~, not_redun_ineq, not_redun_eq] = polyhedron.lcon2vert(A_ineq_both, b_ineq_both, A_eq_both, b_eq_both);
      % not_redun_ineq
      % not_redun_eq

      result = HalfspaceRepresentation(A_ineq_both, b_ineq_both, A_eq_both, b_eq_both);
    end % End of function


    %% Overload "==" operator
    function is_equal = eq(left, right)
      left_vertices = left.getBoxBoundedVertices();
      right_vertices = right.getBoxBoundedVertices();
      left_subset_right = all(ismember( left_vertices', right_vertices', 'rows'));
      right_subset_left = all(ismember(right_vertices',  left_vertices', 'rows'));

      % % TODO: Use a tolerance.
      % pwintz.arrays.findColumnsIn(left_vertices, right_vertices, tolerance=1e-8);
      % pwintz.arrays.findColumnsIn(right_vertices, left_vertices, tolerance=1e-8);

      is_equal = left_subset_right && right_subset_left;
    end % End of function
    
    % Overload "~=" operator
    function is_not_equal = ne(left, right)
      is_not_equal = ~eq(left, right);
    end % End of function

    % ╭───────────────────────────────────────────────────────────────╮
    % │  ╭─────────────────────────────────────────────────────────╮  │
    % │  │             Overload String Representations             │  │
    % │  ╰─────────────────────────────────────────────────────────╯  │
    % ╰───────────────────────────────────────────────────────────────╯
    
    function string_rep = char(this)
      % This method defines the string representation that is 
      % inserted when using "%s" in a string format, such as 
      % fprintf("%s", polyhedron) or sprintf("%s", polyhedron); 
      % It does not change the "disp" output. 
      arguments(Output)
        string_rep char; % Ensure output is cast to char, even if you create a string.
      end

      string_rep = pwintz.strings.format("HalfspaceRepresentation with ambient_dimension=%d, n_inequality_constraints=%d, n_equality_constraints=%d:\n\tA_ineq=%.2g\n\tb_ineq=%.2g\n\tA_eq=%.2g\n\tb_eq=%.2g", this.ambient_dimension, this.n_inequality_constraints, this.n_equality_constraints, this.A_ineq, this.b_ineq, this.A_eq, this.b_eq);
    end

    function disp(this)
      fprintf('%s\n', this);
    end % End of function

    
    function plot(this)
      
      for i_ineq = 1:this.n_inequality_constraints
        % A_ineq * x \leq b_ineq
        a = this.A_ineq(i_ineq, :);
        b = this.b_ineq(i_ineq);
        x0 = b * a' / norm(a);
        a_perp = [0, 1; -1, 0] * a';
        p = [x0 + 3* a_perp,  x0 - 3* a_perp];
        pwintz.plots.plotVector2(x0, 0.5*a);
        plot(x0(1), x0(2), "*");
        plot(p(1, :), p(2, :), "-");
      end


    end

  end % End of methods block


  methods(Access = {?matlab.unittest.TestCase}) % Private block, accessible by unit tests.
    
    function [vertices, is_vertex_on_boundary_of_box] = getBoxBoundedVertices(this)
      % Find the vertices of the intersection of this polyhedron intersected with the bounding box.
      
      % `  1 * x1 \leq 1e6
      % ` -1 * x1 \leq 1e6 
      A_with_bnd_box = [this.A_ineq; this.A_bounding_box];
      b_with_bnd_box = [this.b_ineq; this.b_bounding_box];
      % A_with_bnd_box = [this.A_bounding_box];
      % b_with_bnd_box = [this.b_bounding_box];
    
      try
        tolerance = 1e-6;
        vertices = polyhedron.lcon2vert(A_with_bnd_box, b_with_bnd_box, this.A_eq, this.b_eq, tolerance, verbose=false);
      catch exception
        % unbounded_error_message = "Non-bounding constraints detected. (Consider box constraints on variables.)";
        % if startsWith(exception.message, unbounded_error_message)
        %   is_bounded = false;
        % else
        msg = sprintf("getBoxBoundedVertices(): Unexpected error when finding the vertices of the polyhedron %s.", this);
        exception = exception.addCause(MException("HalfspaceRepresentation:isBounded", msg));
        rethrow(exception);
        % end
      end % End of try-catch block

      n_verts = size(vertices, 2);
      if n_verts >= this.ambient_dimension + 1 % ? How many vertices do we need for convhull to work?
        try
          % Use convhull to ensure that in 2D the vertices are sorted in a counter-clockwise direction.
          cvx_hull_ndxs =  convhull(vertices);
        catch exception
          pwintz.strings.format("convhull failed to compute the convex hull of %d vertices in %dD. vertices = %s", n_verts, this.ambient_dimension, size(vertices, 1))
          rethrow(exception);
        end % End of try-catch block
        % The output "cvx_hull_ndxs" from convhull has the same index at the start and end, so we leave off the end point.
        vertices = vertices(cvx_hull_ndxs(1:end-1), :);
      end 
      
      % Try to fix the twisting of sets. 
%       n_verts = size(vertices, 1);
%       if n_verts >= this.dimension + 1
%       convhulln(vertices)
%       end
      
      % Transpose "vertices" so that each column contains a vertex...
      vertices = transpose(vertices);

      is_vertex_entry_near_boundary = abs(vertices) >= this.bounding_box_half_width * 0.999;

      % If any of the entries of a vertex is at the boundary of the box (within some tolerance), then the vertex is considered on the boundary of the box.
      is_vertex_on_boundary_of_box = any(is_vertex_entry_near_boundary, 1);
       
    end
  end % End private methods block

end % End of class


function mustHaveOnlyNonzeroRows(array)
  % Check that none of the rows are zero. 
  zero_rows = pwintz.arrays.findRowIn(zeros(1, size(array, 2)), array, tolerance=1e-12);
  if ~isempty(zero_rows)
    error("HalfspaceRepresentation:mustHaveOnlyNonzeroRows", "%s=%s must not have a zero row.", mat2str(array));
  end
end % end function

