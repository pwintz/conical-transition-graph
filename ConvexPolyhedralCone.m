classdef ConvexPolyhedralCone < ConvexPolyhedron
	% Tests for this class are located in TestConvexPolyhedralCone and can be run using ConvexPolyhedralCone.test().
	% ` runtests TestConvexPolyhedralCone

	properties(SetAccess = immutable)
		% Define instance constants.
		% dimension (1, 1) int32; % Dimension of ambient space containing the cone.
		% rays (:, :) double; % Non-zero vectors that define the cone as cone(v1, v2, ..., vn);
		% vertex (:, 1); % The origin of the cone.

		% A vector that is strictly inside the cone
		% middle_vector (:, :) double {};
		rays (:, :) double; % Array of rays, given as columns that define the cone.
		n_rays (1, 1) {mustBeInteger,mustBeNonnegative};
	end

	methods(Static)
		function cone = fromRays(rays)
			halfspace_representation = HalfspaceRepresentation.fromConicalHull(rays);
			cone = ConvexPolyhedralCone(halfspace_representation, rays);
		end % End of function

		function empty_array = empty()
			empty_array = ConvexPolyhedralCone();
		end % End of function
	end

	methods
		% Constructor
		function this = ConvexPolyhedralCone(halfspace_representation, rays)
			arguments(Input)
				% Spanning vectors as an array, with each vector a column.
				halfspace_representation (1, 1) HalfspaceRepresentation {mustBeNonempty} = HalfspaceRepresentation(); 
				rays (:, :) double = [];
			end % End of Input arguments block
			this = this@ConvexPolyhedron(halfspace_representation);

			% dimension = halfspace_representation.ambient_dimension;
			% n_rays = size(rays, 2);
			% origin = zeros(this.dimension, 1);

			origin = zeros(this.dimension, 1);
			assert(~pwintz.arrays.isColumnIn(origin, rays, tolerance=1e-9), "The origin should not be in the list of rays.");
			% origin_ndx = pwintz.arrays.findColumnIn(origin, rays, tolerance=1e-10);
			% n_rays = size(rays, 2);
			% rays = rays(:, setdiff(1:n_rays, origin_ndx));
			% assert(isempty(origin_ndx))

			if isempty(rays)
				% disp("Constructing rays using halfspace halfspace_representation.getConeRays()");
				rays = halfspace_representation.getConeRays();
			end
			this.n_rays = size(rays, 2);
			this.rays = pwintz.arrays.normalizeColumns(rays);

			% TODO: Check that the given rays match the given halfspace representation.

			% Check that the polyhdron contains the rays.
			halfspace_representation.containsPoints(rays);

			return

			% ╭──────────────────────────────────────────────────────────────────────────────╮
			% │             Check that vectors are not multiples of each other               │
			% ╰──────────────────────────────────────────────────────────────────────────────╯
			% Normalize the spanning vectors. This is safe because we have asserted that it does not include zero vector.
			% rays = rays ./ vecnorm(rays);
			% rays = pwintz.arrays.uniqueColumns(rays);
			% assert(n_rays == size(rays, 2), "There were redundant vectors in rays");

			% ╭───────────────────────────────────────────────────────────────────────╮
			% │             Check that spanning vectors are not redundant             │
			% ╰───────────────────────────────────────────────────────────────────────╯
			
			% if n_rays >= 3
			%   X = [rays'; origin'];
			%   indices = convhull(X, 'simplify', true);
			% 
			%   % Check that all of the rays are on the boundary of the cone.
			%   assert(numel(indices) == n_rays + 1, "One or more of the given rays were redundant or inside the cone.");
			% elseif n_rays == 2
			%   pwintz.assertions.assertNotEqual(n_rays(:, 1), -n_rays(:, 2), reason="The two rays are antiparallel, producing a full line. We do not support this yet.");
			% end

			% assert(dimension == 2, "Only 2D is implemented. Instead dimension=%d.", dimension);
			
			% TODO: Check that none of the rays are zero. 
			% TODO: Normalize rays. 
			

			
			% ╭────────────────────────────────────────────────────────────────╮
			% │             Construct the halfspace representation             │
			% ╰────────────────────────────────────────────────────────────────╯
			% Use the "Analyze N-dimensional Convex Polyhedra" MATLAB package from File Exchange to convert the polytope into its halfspace-representation.
			% https://www.mathworks.com/matlabcentral/fileexchange/30892-analyze-n-dimensional-convex-polyhedra
%       polytope = [this.rays, origin];
%       [A,b,Aeq,beq] = vert2lcon(polytope');
% 
%       [~, zero_ieq_constraint_ndxs] = ismembertol(0, b, 1e-8, "OutputAllIndices", true);
%       zero_ieq_constraint_ndxs = zero_ieq_constraint_ndxs{1}; % Get result as non-cell array.
%       % zero_ieq_constraint_ndxs = setdiff(1:numel(b), nonzero_ieq_constraint_ndx)
%       A = A(zero_ieq_constraint_ndxs, :);
%       b = b(zero_ieq_constraint_ndxs);
%       this.halfspace_representation = struct("A", A, "b", 0*b, "Aeq", Aeq, "beq", beq);

			% TODO: Check that polyhedron is convex.
			% this.center = mean(rays, 1);
			% 
			% % ⋘────────── Find lexicographically "first" vertex ──────────⋙
			% [~, y_min_ndxs] = min(rays(2,:));
			% if numel(y_min_ndxs) == 1
			%   this.first_vertex_ndx = y_min_ndxs;
			% else
			%   [~, x_min_ndx] = min(rays(1, y_min_ndxs));
			%   assert(numel(x_min_ndx) == 1, "There should be a unique 'first' vertex")
			%   this.first_vertex_ndx = y_min_ndxs(x_min_ndx);
			% end
			% this.first_vertex_ndx 
		end
	end

	methods
		% ╭────────────────────────────────────────────╮
		% │ ╭────────────────────────────────────────╮ │
		% │ │             Set Operations             │ │
		% │ ╰────────────────────────────────────────╯ │
		% ╰────────────────────────────────────────────╯
		function result = minkowskiSum(this, polytope)
			arguments(Input)
				this      %(1, 1)  {pwintz.validators.mustBeScalar}; % ConvexPolyhedralCone
				polytope  %(1, 1)  {pwintz.validators.mustBeScalar}; % Polytope
			end % End of Input arguments block
			
			
			pwintz.strings.format("ConvexPolyhedralCone.minkowskiSum(this=%spolytope=%s)", this, polytope);

			% disp("Using ConvexPolyhedralCone/minkowskiSum");
			% ! We do a simple and inefficient implementation of the Minkowski sum here. 
			% ! It should be fine for small numbers of points, but the complexity is O(m*n) (or worse) instead of O(m + n), where m and n are the numbers of rays in each polygon.
			if class(polytope) ~= "Polytope"
				msg = pwintz.strings.format("Only implemented for sum of a polytope and a convex polyhedral cone. The polytope argument=%s was instead a %s.", polytope, class(polytope));
				error(msg);
			end
			pwintz.assertions.assertEqual(this.dimension, polytope.dimension, leftName="Cone dimension (left argument)", rightName="Polytope dimension (right argument)");

			% ⋘────────── Get numbers of rays ──────────⋙
			n_cone_rays = size(this.rays, 2);
			n_polytope_verts = size(polytope.vertices, 2);
			n_points = n_cone_rays * n_polytope_verts;

			% ⋘────────── Generate all the combinations of sums ──────────⋙
			points = nan(this.dimension, n_points);
			i_point = 1;
			for i_ray = 1:n_cone_rays
				for i_vertex = 1:n_polytope_verts
					% Brute force addition of all vertices and rays.
					points(:, i_point) = this.rays(:, i_ray) + polytope.vertices(:, i_vertex);
					i_point = i_point + 1;
				end
			end

			% Since we are adding a cone, which includes the origin, all of the points in the polytope should also be included. 
			points = [polytope.vertices, points];

			% Construct a halfspace for the restricted set, formed frmo using only the finite representation of the rays. 
			% After, we must remove any half space restrictions defined by planes that don't touch any of the vertices in the original polytope "polytope".
			prelim_halfspace_rep = HalfspaceRepresentation.fromConvexHull(points);

			A_ineq = prelim_halfspace_rep.A_ineq;
			b_ineq = prelim_halfspace_rep.b_ineq;
			A_eq = prelim_halfspace_rep.A_eq;
			b_eq = prelim_halfspace_rep.b_eq;

			% Check that all of the polytope's vertices are in the Minkowski sum.
			assert(all(A_ineq * polytope.vertices - b_ineq <= 1e-8, "all"));

			% ⋘──────── Find spurious inequality constraints that were introduced by the bounding box. ────────⋙
			% Create a matrix where each (i,j) entry is true if the jth vertex of the polytope is on the boudnary of the ith half space.
			% We do this by evaluating 
			%   A * [p1, p2, ... pn] - b      (*)
			% where p1, p2, etc. are the vertices of the polytope.  
			% If the (i, j) entry of (*) is near zero, then the p_j vertex is near the boundary of the ith half space. 
			tolerance_to_be_in_bnd = 1e-6;
			is_polytope_vert_on_halfspace_bnd = abs(A_ineq * polytope.vertices - b_ineq) < tolerance_to_be_in_bnd;

			% Create a column vector such that each ith entry is 1 if none of the vertices in P are in the boudnary of the jth halfspace. That is, if all of the entries in the jth row are zero.
			% is_halfspace_spurious = ~any(is_polytope_vert_on_halfspace_bnd, 2); % Any "in each row".
			is_halfspace_spurious = pwintz.logical.isNoneForEachRow(is_polytope_vert_on_halfspace_bnd);

			% Delete spurious half spaces. 
			A_ineq = A_ineq(~is_halfspace_spurious, :);
			b_ineq = b_ineq(~is_halfspace_spurious);

			halfspace_rep = HalfspaceRepresentation(A_ineq, b_ineq, A_eq, b_eq);

			% ⋘────────── Set result to convex hull ──────────⋙
			result = ConvexPolyhedron(halfspace_rep);

			% # The following code was an attempt at a more efficient implementation.
			% https://cp-algorithms.com/geometry/minkowski.html
			% left_ndx  = left.first_vertex_ndx;
			% right_ndx = right.first_vertex_ndx;
			% result_n_rays = left.n_rays + right.n_rays;
			% for result_ndx = 1:result_n_rays
			%   sult_n_rays(result_ndx) = left.rays(left_ndx) + right.rays(right_ndx);
			%   if atan2(y, x)
			% end
		end

		function result = linearTransform(this, A)
			arguments(Input)
				this ConvexPolyhedralCone;
				A    (:, :) double;
			end % End of Input arguments block
			% A
			% this.rays
			transformed_rays = A * this.rays;
			result = ConvexPolyhedralCone.fromRays(transformed_rays);
		end

		% ╭─────────────────────────────────────────────────────────────╮
		% │             Override intersection for two cones             │
		% ╰─────────────────────────────────────────────────────────────╯
		function result = intersection(this, other)
			% The intersection of two cones is also a cone, so we 
			if isa(other, "ConvexPolyhedralCone")

				%this_hs  = this.halfspace_representation;
				%other_hs = other.halfspace_representation;
				%this_rays  = this.rays;
				%other_rays = other.rays;
				%% is_this_rays_in_other = this.halfspace_representation.containsPoint
				%is_other_rays_in_this = this_hs.containsPoints(other_rays);
				%is_this_rays_in_other = other_hs.containsPoints(this_rays);
				%this_rays_in_other = this_rays(:, is_this_rays_in_other);
				%other_rays_in_this = other_rays(:, is_other_rays_in_this);
				%combined_rays = [this_rays(:, is_this_rays_in_other), other_rays(:, is_other_rays_in_this)];
				intersection_halfspace_representation = intersection(this.halfspace_representation, other.halfspace_representation);
				result = ConvexPolyhedralCone(intersection_halfspace_representation);
			else
				result = intersection@ConvexPolyhedron(this, other);
			end
		end

		% ╭─────────────────────────────────────╮
		% │             Unit Sphere             │
		% ╰─────────────────────────────────────╯
		function sphere_over_approx_poly = getUnitSphereOverApproximation(this)
		  max_angle = this.getMaxAngle();
		  assert(max_angle < pi, 'max_angle must be less than pi.');
		  % From geometry we have that convex hull of {v1, v2, ..., vp, r*v1, r*v2, ..., r*vp} contains the unit sphere if r = sec(theta_max / 2).
		  r = sec(max_angle / 2); 
			
			% if size(this.rays, 1) == 2
			% 	% Sort the rays
			% 	% vertices = [this.rays(:, 1), r*this.rays, this.rays(:, 2)]
			% 	vertices = [this.rays, r*this.rays];
			% 	sphere_over_approx_poly = Polytope.fromConvexHull(vertices);
			% else
				sphere_over_approx_poly = Polytope.fromConvexHull([this.rays, r*this.rays]);
			% end
		end

		function max_angle = getMaxAngle(this)
		  % Compute the largest angle between rays in this cone.
		  % To find the largest angle, we use the fact that the angle theta between the vertices v_i and v_j (which are unit vectors) is given by cos(theta) = v_i'*v_j. The value of theta is maximized for the choice of (i,j) such that v_i'*v_j is minimized.
		  % if this.dimension ~= 2
		  %   warning('Only tested for 2 dimensions, but it should work.');
		  % end
		  % [~, adjacent_vertices] = this.getVerticesAdjacentToCone(cone_ndx);
		  % Compute V'*V, where V = [v1, v2, ..., vp] is the matrix formed from all the cones boundary rays. 
		  dot_products_between_each_ray = this.rays'*this.rays;
		  min_dot_product = min(dot_products_between_each_ray, [], 'all');
		  max_angle = acos(min_dot_product);
		end

		% ╭──────────────────────────────────────────────────────╮
		% │ ╭──────────────────────────────────────────────────╮ │
		% │ │             Overload Default Methods             │ │
		% │ ╰──────────────────────────────────────────────────╯ │
		% ╰──────────────────────────────────────────────────────╯
		% % ╭───────────────────────────────────────────────────────────────────────╮
		% % │             Make ConvexPolyhedralCone Behave Like a Container             │
		% % ╰───────────────────────────────────────────────────────────────────────╯
		% function result = isempty(this)
		% 	result = isempty(this.rays);
		% end
		% function result = count(this)
		% 	result = size(this.rays, 2);
		% end

		% ╭────────────────────────────────────────────╮
		% │             Overload Operators             │
		% ╰────────────────────────────────────────────╯
		% % Overload "|" operator to compute the intersection of sets.
		% function result = or(left, right)
		% 	result = intersection(left, right);
		% end

		% Overload "+" operator to compute the Minkowski sum of sets.
		function result = plus(left, right)
			result = minkowskiSum(left, right);
		end

		% Overload "*" operator
		function result = mtimes(A, convex_polyhedral_cone)
			result = convex_polyhedral_cone.linearTransform(A);
		end

		% Overload "==" operator
		function is_equal = eq(left, right)
			left_subset_right = all(ismember( left.rays', right.rays', 'rows'));
			right_subset_left = all(ismember(right.rays',  left.rays', 'rows'));
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
			if this.n_rays == 0
				rays_string = sprintf('0 rays (empty).\n');
			elseif this.n_rays == 1
				rays_string = sprintf('1 ray: \n%s\n', formattedDisplayText(this.rays));
			else
				rays_string = sprintf('%d rays: \n%s\n', this.n_rays, formattedDisplayText(this.rays));
			end
			str = sprintf('ConvexPolyhedralCone in %dD with %sand with halfspace representation: \n%s\n', this.dimension, rays_string, this.halfspace_representation);
			% if this.n_rays == 0
			%   str = sprintf('ConvexPolyhedralCone with 0 rays (empty).\n');
			% elseif this.n_rays == 1
			%   str = sprintf('ConvexPolyhedralCone with 1 spanning_vector: \n%s\n', formattedDisplayText(this.rays));
			% else
			%   str = sprintf('ConvexPolyhedralCone with %d rays: \n%s\n', this.n_rays, formattedDisplayText(this.rays));
			% end
		end

		function disp(this)
			if ~isscalar(this)
				fprintf('%dx%d ConvexPolyhedralCone array\n', size(this, 1), size(this, 2))
				return
			end
			fprintf("%s\n", char(this));
		end

		function out = isequaln(left, right)
			out = eq(left, right);
		end % End of function

		% ╭────────────────────────────────────────╮
		% │  ╭──────────────────────────────────╮  │
		% │  │             Plotting             │  │
		% │  ╰──────────────────────────────────╯  │
		% ╰────────────────────────────────────────╯
		function plot(this, varargin)
			color = 0.0 * [1, 1, 1];
			switch this.dimension
				case 2
					pts = [this.rays, 0*this.rays(:, 1)];
					patch(pts(1,:), pts(2, :), color, "FaceAlpha", 0.2, varargin{:});
				case 3
					edges = [
						1, 2, 3;
						1, 2, 4;
						1, 3, 4;
					];
					% points = [0*this.rays(:, 1), this.rays]'
					% plot(pts(1,:), pts(1, :), pts(3,:));
					% tetramesh(edges, points, color, "FaceAlpha", 0.2, varargin{:});
					this.rays
					origin = zeros(this.dimension, 1);
					vecnorm(this.rays)
					patch(this.rays(1,:), this.rays(2, :), this.rays(3, :), color, "FaceAlpha", 0.2, varargin{:});
					for i = 1:(this.n_rays-1)
						points = [this.rays(:, i), this.rays(:, i+1), origin];
						patch(points(1,:), points(2, :), points(3, :), color, "FaceAlpha", 0.2, varargin{:});
					end
				otherwise
					error("Unexpected case: %s.", this.dimension);
			end
			
			
		end % End of function
	end% End of methods block

	% methods(Access=private, Static)
	% 	function rays = raysFromOperand(operand)
	% 		% Given an operand that is either a ConvexPolyhedralCone or a 2xN array, return the corresponding set of rays. For the case where "operand" is an array, the "rays" are the "operand" array.
	% 		
	% 		if isa(operand, "ConvexPolyhedralCone")
	% 			rays = operand.rays;
	% 		elseif isnumeric(operand) && size(operand, 1) == 2
	% 			rays = operand;
	% 		else
	% 			error("ConvexPolyhedralCone:invalidOperand","operand is invalid type (%s) or size (%s). Expected a Convex Polyhedron or a 2xN array.", class(operand), mat2str(size(operand)));
	% 		end
	% 	end
	% end % End of private static methods block

	% methods(Access=protected)
	%   function groups = getPropertyGroups(obj)
	%       if isscalar(obj)
	%         % Specify the values to be displayed for properties
	%         propList = struct('rays', newline + formattedDisplayText(obj.rays));%) mat2str(obj.rays, 3));
	%         groups = matlab.mixin.util.PropertyGroup(propList);
	%       else
	%         % Nonscalar case: call superclass method
	%         groups = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
	%       end
	%   end
	% end

	% TODO: Try to make nicer formatting for the strings.
	methods%(Access=private)
		% function str = raysToString(this)
		%   string_array = strings(2, this.n_rays)
		% end
	end
end


