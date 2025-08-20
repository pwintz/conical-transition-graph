classdef ConicalPartition < handle
  %%% ConicalPartition - Create a class for defining a partition of R^n into cones. 
  %%% The cones are defined usings a "watermelon slicing" method using equal angles for all cones.
  %%% To keep things simple, we are starting with 2D. 
  % 
  % Unit tests can be run via ConicalPartition.test().

  properties(Constant, GetAccess = {?matlab.unittest.TestCase}) % Define class constants.
    % Define the smallest angle between vertex angles before two vertices are merged. 
    ANGLE_TOL (1, 1) double = 1e-5;
  end % End of class constants block

  properties(SetAccess = immutable)
    % Define instance constants.
    A (2,2) double; % Matrix of "\dot x = A x".
    dimension (1,1);  % The dimension of the system considered
    
    derivative_vertex_angles (1, :) double; 
    state_vertex_angles      (1, :) double;     
    derivative_vertices      (2, :) double;     
    state_vertices           (2, :) double; 
  end

  properties(SetAccess = immutable, GetAccess=private)
    % n_slices (1,1); % Number of divisions per dimension  
    % derivative_space_angle; % The angle swept by each cone in derivative space along the rotation in each dimension.
    % state_vertex_angles; % The angle of each slice in state space.
    vertex_multi_index_size (1,1);
  end

  properties(Dependent)
    % Define dependent variables.
    n_vertices (1,1); % The numbers of state vertices and derivative vertices are equal.
    vertex_indices (1, :) int32;
    region_indices (1, :) int32;
  end
  methods
    function value = get.n_vertices(this)
      value = size(this.state_vertices, 2);
    end
    function value = get.vertex_indices(this)
      value = 1:this.n_vertices;
    end
    function value = get.region_indices(this)
      value = 1:this.n_vertices;
    end
  end

  methods(Static) % Define convenience functions for running tests.
    function test(test_function_name)
      test_name = "TestConicalPartition";
      if nargin() == 1 % If a test function name is given, append it so only that test is run.
        test_name = test_name + "/" + test_function_name;
      end
      results = runtests(test_name);
      fprintf("%d Passed, %d Failed, %d Incomplete.\n", sum([results.Passed]), sum([results.Failed]), sum([results.Incomplete]));
    end % End of function
  end % End static methods block


  methods(Static)
    function conical_partition = fromNumberOfUniformDerivativeSlices(A, n_slices)
      arguments(Input)
        A (2, 2) double;
        n_slices (1,1) double;
      end % End of Input arguments block
      
      arguments(Output)
        conical_partition ConicalPartition;
      end % End of Output arguments block
      
      derivative_vertex_angles = linspace(0, 2*pi, n_slices + 1);
      derivative_vertex_angles = derivative_vertex_angles(1:end-1);
      conical_partition = ConicalPartition(A, derivative_vertex_angles);
    end

    function conical_partition = fromAngles(A, derivative_vertex_angles, extra_state_vertex_angles)
      % p = inputParser;
      % addRequired(p, "A", @(A) isnumeric(A) && ismatrix(A) && issquare(A)); % Positional argument.
      % addParameter(p, "nSlices", 8, @(x) isreal(x) && rem(x,1)==0 && x > 1); % Name-value parameter.;
      % parse(p, varargin{:});

      if isempty(extra_state_vertex_angles)
        combined_derivative_vertex_angles = derivative_vertex_angles;
      else
        extra_state_vertices = pwintz.math.angle2UnitVector(extra_state_vertex_angles);
        extra_derivative_vertex_angles = pwintz.math.atan2(A * extra_state_vertices);

        % Combine the given derivative vertex angles with the ones coming from the extra state vertex angles. We don't need to sort them here because that is done in the ConicalPartition constructor.
        combined_derivative_vertex_angles = [derivative_vertex_angles, extra_derivative_vertex_angles];
      end

      % vertices_table = struct2table(struct(...
      %     "state_vertex",      state_vertices, ...
      %     "derivative_vertex", derivative_vertices, ...
      %     "state_angle",       state_vertex_angles, ...
      %     "derivative_angle",  derivative_vertex_angles  ...
      % ))
      conical_partition = ConicalPartition(A, combined_derivative_vertex_angles);
    end % End of function
  end % End static methods block

  methods(Access = private)

    % Constructor
    function this = ConicalPartition(A, derivative_vertex_angles)
      arguments(Input)
        A (2,2) double;
        derivative_vertex_angles (1, :) double;
      end % End of Input arguments block
      
      % Parse inputs from varargin argument.
      % p = inputParser;
      % addRequired(p, "A", @(x) isnumeric(x) && ismatrix(x)); % Positional argument.
      % addParameter(p, "dimension", 2, @(x) x == 2); % Name-value parameter.
      % addParameter(p, "nSlices", 8, @(x) isreal(x) && rem(x,1)==0 && x > 1); % Name-value parameter.;
      % addParameter(p, "stateAngles", []); % Name-value parameter.;
      
      % parse(p, varargin{:});
      % args = p.Results;
      % this.A = args.A;
      % this.dimension    = args.dimension;
      % this.n_slices      = args.nSlices;

      % ╭───────────────────────────────────────────────╮
      % │             Check Input Arguments             │
      % ╰───────────────────────────────────────────────╯
      % ⋘────────── Check "A" ──────────⋙
      assert(isnumeric(A));
      assert(isreal(A));
      assert(ismatrix(A));
      assert(size(A,1) == size(A,2), "The matrix A must be square");
      assert(cond(A) < 1e9, "The matrix A must be well-conditioned but had a condition number of %8.2g", cond(this.A))
      
      assert(size(A, 1) == 2, "Only 2D implemented, currently.");

      % ⋘────────── Check "derivative_vertex_angles" ──────────⋙
      

      % ╭────────────────────────────────────────────────╮
      % │             Assign Property Values             │
      % ╰────────────────────────────────────────────────╯
      this.A = A;
      this.dimension = size(A, 1); % Equal to size(A, 2) by assertion.
      % this.n_slices      = args.nSlices;
      this.vertex_multi_index_size = 1; % For 2D Case. Must be assigned before the first use of get...Vertex()

      % ⋘────────── Sort and Filter Duplicate Derivative Angles ──────────⋙
      derivative_vertex_angles = mod(derivative_vertex_angles, 2*pi); % Map to [0, 2*pi)
      derivative_vertex_angles = uniquetol(derivative_vertex_angles, ConicalPartition.ANGLE_TOL);
      derivative_vertex_angles = sort(derivative_vertex_angles);
% 
%       % ⋘──── Compute Derivative Vertex Angles of Extra State Vertex Angles ────⋙
%       for i = 1:numel(extra_state_vertex_angles)
%         extra_state_vertices(:,i)           = [cos(extra_state_vertex_angles(i)); sin(extra_state_vertex_angles(i))];
%         extra_derivative_vertex_angles(:,i) = pwintz.math.atan2(A * extra_state_vertices(:,i));
%       end
% 

      % ⋘──── Compute Property Values Derived from Given Values ────⋙
      % n_slices = numel(derivative_vertex_angles);
      % derivative_vertices  = nan([this.dimension, n_slices]);
      % state_vertices       = nan([this.dimension, n_slices]);


      derivative_vertices = pwintz.math.angle2UnitVector(derivative_vertex_angles);
      state_vertices      = A \ derivative_vertices;
      state_vertex_angles = pwintz.math.atan2(state_vertices);

      % Remap from [-pi, pi] to [0, 2pi].
      state_vertex_angles = mod(state_vertex_angles, 2*pi); 

      this.derivative_vertices      = derivative_vertices;
      this.state_vertices           = state_vertices;
      this.derivative_vertex_angles = derivative_vertex_angles;
      this.state_vertex_angles      = state_vertex_angles;

%       for i = 1:n_slices
%         % derivative_vertices(:,i) = [cos(derivative_vertex_angles(i)); sin(derivative_vertex_angles(i))];
%         state_vertices(:,i)      = A \ derivative_vertices(:,i);
%         % ⋘────────── Compute state vertex angle ──────────⋙
%         state_angle = pwintz.math.atan2(state_vertices(:,i));
%         state_angle = mod(state_angle, 2*pi); % Remap from [-pi, pi] to [0, 2pi].
% 
%         state_vertex_angles(i) = state_angle;
%       end





      % ! Switching to a table as the data type to store vertices might be easier.
      % derivative_vertices = struct2table(struct(...
      %     "x", value1, ...
      %     "y", value1, ...
      %     "angle", value2  ...
      %   ) ...
      % );

      % ⋘────────── Check Properties ──────────⋙

      ConicalPartition.validateAngles(this.derivative_vertex_angles, this.state_vertex_angles);
      % It is important that derivative_vertex_angles is sorted because we use that property in getDerivativeRegionIndexContainingPoint.
      % assert(all(this.derivative_vertex_angles == sort(this.derivative_vertex_angles)), "derivative_vertex_angles must be sorted.");
      % assert(...
      %   all(diff([this.state_vertex_angles, 2*pi]) <= pi), ...
      %   "All of the state vertex angles must be less than 180 degress so that the cones are convex.\n\tstate_vertex_angles=%s", mat2str(this.state_vertex_angles)...
      % );
      % assert(...
      %   all(diff([this.derivative_vertex_angles, 2*pi]) <= pi), ...
      %   "All of the state vertex angles must be less than 180 degress so that the cones are convex."...
      % );
    end
  end % End private methods

  methods % Public methods
    % ╭──────────────────────────────────╮
    % │             Indexing             │
    % ╰──────────────────────────────────╯
    function multi_indexArray = getVertexMultiIndices(this)
      % vertexMultiIndices - Get an wide array where each column is the multi-index of a one vertex in this conical partition.
      % The result can be used like 
      % *
      % *  for i = conical_partition.getVertexMultiIndices()
      % *    v_i = conical_partition.getStateSpaceVertex(i)
      % *  end
      % *
      assert(this.dimension == 2, "Only 2D implemented.");

      % For the 2D case, each multi_indexing only has one entry.
      multi_indexArray = 1:this.n_vertices; %  <- For 2D.
      % multi_indexArray = 1:this.n_slices;
    end
    
    function v = getStateSpaceVertex(this, multi_index)
      % getVertex - Given the multi-index of the vertex, passed as the arguments, returns the vertex as a column vector.
      vertex_index = this.normalizeIndex(multi_index);
      v = this.getDerivativeSpaceVertex(vertex_index);
      v = this.A \ v;
      v = v / norm(v);
    end

    function v = getDerivativeSpaceVertex(this, multi_index)
      % getVertex - Given the multi-index of the vertex, passed as the arguments, returns the vertex as a column vector.
      vertex_index = this.normalizeIndex(multi_index);
      vertex_angle = this.derivative_vertex_angles(vertex_index);
      v = [cos(vertex_angle); sin(vertex_angle)];

    end

    function theta = getStateSpaceVertexAngle(this, multi_index)
      % getVertex - Given the multi-index of the vertex, passed as the arguments, returns the vertex as a column vector.
      vertex_index = this.normalizeIndex(multi_index);
      theta = this.state_vertex_angles(vertex_index);
    end
    
    function C = getStateCone(this, multi_index)
      % Given an index, get the state space cone for a given multi-index, given as a polyhedron (the cone truncated to the unit sphere), which is encoded as a row of vertices that defines the polyhedron.
      vertex_index = this.normalizeIndex(multi_index);
      v1 = this.getStateSpaceVertex(vertex_index);
      v2 = this.getStateSpaceVertex(vertex_index + 1);
      origin = [0; 0];
      C = ConvexPolyhedron.fromConvexHull([v1, origin, v2]);
    end

    function D = getDerivativeCone(this, multi_index)
      % Given an index, get the state space cone for a given multi-index, given as a polyhedron (the cone truncated to the unit sphere), which is encoded as a row of vertices that defines the polyhedron.
      vertex_index = this.normalizeIndex(multi_index);
      v1 = this.getDerivativeSpaceVertex(vertex_index);
      v2 = this.getDerivativeSpaceVertex(vertex_index+1);
      origin = [0; 0];
      D = ConvexPolyhedron.fromConvexHull([v1, origin, v2]);
    end

    % ╭──────────────────────────────────────────────────────────────────╮
    % │             Translating between vertices and regions             │
    % ╰──────────────────────────────────────────────────────────────────╯
    
    function v_neighbors_ndxs = getNeighborVertexIndices(this, multi_index)
      % getNeighborVertexIndices - Get an wide array containing where each column is the multi-index of a vertex that is adjacent to the given multi-index.
      % The result can be used like 
      % *                                                
      % * for i = conical_partition.getVertexMultiIndices() 
      % *   v_i = conical_partition.getStateSpaceVertex(i)         
      % *   for j = conical_partition.getNeighborVertexIndices(i) 
      % *     v_j = conical_partition.getStateSpaceVertex(i)  
      % *   end                                                   
      % * end                                                   
      % *                                                

      neighborKernel = [-1, +1]; % For 2D case

      v_neighbors_ndxs = this.normalizeIndex(multi_index + neighborKernel);
    end

     function adjacent_regions_ndxs = getAdjacentRegionsIndices(this, v_ndx)
      % getAdjacentRegionsIndices - Get a list of multi-indices for regions with the given vertex (defined by v_ndx) in their boundary.
      % The result is wide array where each column is the multi-index of a region.
      neighborRegionsNdxOffsets = [-1, 0]; % For 2D case

      adjacent_regions_ndxs = this.normalizeIndex(v_ndx + neighborRegionsNdxOffsets);
    end

    function conjoining_regions_ndxs = getConjoiningRegionsIndices(this, v1_ndx, v2_ndx)
      % getConjoiningRegionsIndices - Get a list of regions that are adjacent to two given vertices.
      % The result is wide array where each column is the multi-index of a region.
      v1_adjacent_regions_ndxs = this.getAdjacentRegionsIndices(v1_ndx);
      v2_adjacent_regions_ndxs = this.getAdjacentRegionsIndices(v2_ndx);

      conjoining_regions_ndxs = intersect(v1_adjacent_regions_ndxs, v2_adjacent_regions_ndxs);
    end

    function bnd_v_ndxs = getBoundaryVertices(this, region_ndx)
      assert(this.dimension == 2, 'Only implemented for 2 dimensions');
      bnd_vertex_ndx_offsets = [0, 1];
      bnd_v_ndxs = this.normalizeIndex(region_ndx + bnd_vertex_ndx_offsets);
    end % End of function
    
    function region_ndxs = getRegionsBetweenVerticesFromAngles(this, start_angle, end_angle)
      assert(all(size(start_angle) == [1, 1]), ...
        "start_angle must have shape %s but instead had shape %s", mat2str(size([1, 1])), mat2str(size(start_angle)));
      assert(all(size(end_angle) == [1, 1]), ...
        "end_angle must have shape %s but instead had shape %s", mat2str(size([1, 1])), mat2str(size(end_angle)));

      tol = 1e-10;
      start_vertex_ndx  = find(pwintz.math.angleDistance(this.state_vertex_angles, start_angle) < tol);
      end_vertex_ndx    = find(pwintz.math.angleDistance(this.state_vertex_angles,   end_angle) < tol);
      if numel(start_vertex_ndx) > 1
        error("Found multiple elements in this.state_vertex_angles=%s that match start_angle=%f, namely at indices=%s", mat2str(this.state_vertex_angles), start_angle, mat2str(start_vertex_ndx))
      end
      if numel(end_vertex_ndx) > 1
        error("Found multiple elements in this.state_vertex_angles=%s that match end_angle=%f, namely at indices=%s", mat2str(this.state_vertex_angles), end_angle, mat2str(end_vertex_ndx))
      end

      assert(all(size(start_vertex_ndx) == [1, 1]), ...
        "start_vertex_ndx must be 1x1 but instead had shape %s", mat2str(size(start_vertex_ndx)));
      assert(all(size(end_vertex_ndx) == [1, 1]), ...
        "end_vertex_ndx must be 1x1 but instead had shape %s", mat2str(size(end_vertex_ndx)));
      assert(start_vertex_ndx ~= end_vertex_ndx);
      
      if start_vertex_ndx < end_vertex_ndx
        region_ndxs = start_vertex_ndx:(end_vertex_ndx - 1);
      else
        % Handle the case where the indices loop past 0.
        region_ndxs = [start_vertex_ndx:this.n_vertices, 1:(end_vertex_ndx - 1)];
      end
    end

    % ╭──────────────────────────────────────────╮
    % │             Set Computations             │
    % ╰──────────────────────────────────────────╯

    function region_ndxs = getStateRegionIndexContainingPoint(this, x)
      % Find all of the regions that contain x. 
      assert(this.dimension == 2, 'Only implemented for 2 dimensions');
      assert(iscolumn(x) && numel(x) == 2);

      x_dot = this.A * x;
      region_ndxs = this.getDerivativeRegionIndexContainingPoint(x_dot);
    end % End of function

    function region_ndxs = getDerivativeRegionIndexContainingPoint(this, x_dot)
      % Find all of the regions that contain x_dot. 
      assert(this.dimension == 2, 'Only implemented for 2 dimensions');
      assert(iscolumn(x_dot) && numel(x_dot) == 2);

      % If x_dot is at the origin, then it is in all of the regions.
      if  norm(x_dot) == 0
        region_ndxs = this.region_indices;
        return
      end % End of if " norm(x_dot) == 0" block

      x_dot_angle = pwintz.math.atan2(x_dot);
      x_dot_angle = mod(x_dot_angle, 2*pi); % Convert range to [0, 2pi] instead of [-pi, pi].
      assert(x_dot_angle >= 0);

      % The values derivative_vertex_angles are sorted, by definition, so if we take the last index where the angle is greater than the vertex angle, then we'll have that it is the region where x_dot is between that vertex and the next vertex.
      region_ndxs = find(x_dot_angle >= this.derivative_vertex_angles, 1, 'last');
    end % End of function

    function [region_ndxs, intersections] = getStateRegionIntersections(this, poly)
      arguments(Input)
        this ConicalPartition;
        poly ConvexPolyhedron; % <- We should make a new data type for polys.
      end
      arguments(Output)
        region_ndxs   (1,:); 
        intersections (1,:) cell;
      end

      does_region_intersect     = zeros(size(this.region_indices));
      intersections = cell(size(this.region_indices));

      for region_ndx = this.region_indices
        region = this.getStateCone(region_ndx);
        intersection  = 1e5*region | poly;
        intersections{region_ndx} = intersection;
        does_region_intersect(region_ndx) = ~isempty(intersection);
      end

      % Convert from a logical array to indices.
      % TODO: when region indices are multi-indices, we might need to convert the output of find, which gives linear indices, using ind2sub, or something like that.
      region_ndxs = find(does_region_intersect);
    end % End of function

    % % ╭──────────────────────────────────────╮
    % % │             Modification             │
    % % ╰──────────────────────────────────────╯
    % function new_conical_partition = refineViaNewStateRegion(this, P)
    %   
    % end % End of function

    % ╭────────────────────────────────────────╮
    % │             Reachable Sets             │
    % ╰────────────────────────────────────────╯
    function R = computeReachableSet(this, region_multi_index, P0)
      % Compute the region that is reachable from P0 when flowing according to \dot x \in AD, where D is the derivative cone for the given multi_index.
      D = conical_partition.getDerivativeCone(region_multi_index);
      % The "derivative cone" is represented using vertices on the unit sphere, so we scale it up so that we get a set the sufficiently represents "P0 + AD", at least up to a large radius.
      R = addDerivativeCone(P0, 1000*D);
      C = conical_partition.getStateCone(region_multi_index);
      R_in_C = intersection(R, 1000*C);
    end

  end  % End of methods


  methods(Access = {?matlab.unittest.TestCase}) % Define private class methods.

    function out = normalizeIndex(this, multi_index)
      % normalizeIndex - Given a multi-index, normalize its entries by applying the modulo operation so that each entry is in {1, 2, ..., n_slices}.

      % Check that the multi-index matches the dimension of the system.
      assert(~isempty(multi_index), "multi_index must be nonempty.");
      assert( ...
        size(multi_index, 1) == this.vertex_multi_index_size, ...
        "The size(multi_index, 1) = %d must be vertex_multi_index_size = %d.", size(multi_index, 1), this.vertex_multi_index_size ...
      )

      % Compute the index of the vertex from the given argument, modulo the number of slices. Since MATLAB uses 1-based indexing, we need to subtract 1 before computing the modulo and then add it afterward.
      assert(this.dimension == 2, "Only 2D implemented.");
      out = mod(multi_index - 1, this.n_vertices) + 1;
    end
  end% End of private methods

  methods(Static, Access = {?matlab.unittest.TestCase}) % Some static functions for internal use.
    
      function validateAngles(derivative_vertex_angles, state_vertex_angles)

        % ╭────────────────────────────────────────────────────────╮
        % │             Check derivative_vertex_angles             │
        % ╰────────────────────────────────────────────────────────╯
        % ⋘────────── Check range ──────────⋙
        if any(derivative_vertex_angles < 0 | derivative_vertex_angles >= 2*pi)
          out_of_range_ndxs = find(derivative_vertex_angles < 0 | derivative_vertex_angles >= 2*pi);
          error(...
            "ConicalPartition:DerivativeAnglesOutOfRange", ...
            'The values of derivative_vertex_angles must in [0, 2*pi). The values %s are out of range at index(s) %s', ...
            mat2str(out_of_range_ndxs), ...
            mat2str(find(out_of_range_ndxs))...
          );
        end

        % ⋘────────── Check strictly increasing ──────────⋙
        if any(diff(derivative_vertex_angles) <= 0)
          error(...
            "ConicalPartition:DerivativeAnglesNonincreasing", ...
            'The values of derivative_vertex_angles must be increasing. The angle decrease at index(s) %s', ...
            mat2str(find(diff(derivative_vertex_angles) > 0))...
          );
        end
        
        % ⋘────────── Check step size ──────────⋙
        % Since the angles "loop around" at the end, we append 2*pi
        derivative_angle_steps = diff([derivative_vertex_angles, derivative_vertex_angles(1) + 2*pi]);
        if any(derivative_angle_steps > pi)
          big_step_ndxs = find(derivative_angle_steps > pi);
          error(...
            "ConicalPartition:DerivativeAngleStepTooBig", ...
            'The difference between consecutive angles in derivative_vertex_angles was %s > pi at index(s). ', ...
            mat2str(derivative_angle_steps(big_step_ndxs)), ...
            mat2str(big_step_ndxs) ...
          );
        end

        % ╭───────────────────────────────────────────────────╮
        % │             Check state_vertex_angles             │
        % ╰───────────────────────────────────────────────────╯
        if nargin() == 2 % The second "state_angle_steps" can be omitted for testing purposes.
          % ⋘────────── Check length matches derivative_vertex_angles ──────────⋙
          assert(all(size(derivative_vertex_angles) == size(state_vertex_angles)));

          derivative_vertex_angles
          state_vertex_angles
          % ⋘────────── Check that there are no duplications ──────────⋙
          [unique_angles, ndxs_of_unique_angles] = uniquetol(state_vertex_angles, ConicalPartition.ANGLE_TOL, 'OutputAllIndices', true);
          if(numel(unique_angles) < numel(state_vertex_angles)) % If there were duplicates:
            % Generate an informative error message...
            is_entry_duplicate = cellfun(@(entry) numel(entry) > 1, ndxs_of_unique_angles');
            duplicate_entries_str = cellfun(@(ndxs) sprintf("%f is duplicated at indices %s", state_vertex_angles(ndxs(1)), mat2str(ndxs)), ndxs_of_unique_angles(is_entry_duplicate))

            % Raise error.
            error("There was a repeated entry in 'state_vertex_angles' (up to a small numerical tolerance). The list had %d items but only %d were unique. List of duplicates: \n\t%s.", numel(state_vertex_angles), numel(ndxs_of_unique_angles), join(duplicate_entries_str, "\n\t"));
          end
          
          % ! We don't need to test for state angles > 180 degrees because this is implied the linear relationship with the derivative angles and the fact that those angles are less than 180 degres.
          % state_angle_steps = pwintz.math.angleDistance(state_vertex_angles, circshift(state_vertex_angles, -1));
          % if any(state_angle_steps > pi)
          %   big_step_ndxs = find(state_angle_steps > pi);
          %   error(...
          %     "ConicalPartition:StateAnglesTooWide", ...
          %     'The difference between consecutive angles in state_vertex_angles was %s > pi at index(s). ', ...
          %     mat2str(state_angle_steps(big_step_ndxs)), ...
          %     mat2str(big_step_ndxs) ...
          %   );
          % end
        end
        
      end % End of function

  end % End static methods block

end%  End of class

% ╭───────────────────────────────────────────────╮
% │  ╭─────────────────────────────────────────╮  │
% │  │             Local functions             │  │
% │  ╰─────────────────────────────────────────╯  │
% ╰───────────────────────────────────────────────╯
function result = addDerivativeCone(P, D)
  % Given a convex polygon P, add the derivative cone D, which consists of three points [v1, 0, v2]. 
  % The addition is in the sense of the Minkowski sum.
  v1 = D(:, 1);
  v2 = D(:, 3);
  result = [P, P + v1, P + v2];
  d_triangulated = delaunayTriangulation(result');
  convex_hull_indices = convexHull(d_triangulated);
  result = d_triangulated.Points(convex_hull_indices,:)';
  % Use geom3d trim extra points
  % result = convexHull(result);
end

function intersection_vertices = intersection(vertices1, vertices2)
  uniqueCols = @(array) unique(array', 'rows', 'stable')';
  vertices1 = uniqueCols(vertices1);
  vertices2 = uniqueCols(vertices2);
  poly1 = polyshape(vertices1(1,:), vertices1(2,:));
  poly2 = polyshape(vertices2(1,:), vertices2(2,:));
  poly_intersection = intersect(poly1, poly2);
  intersection_vertices = poly_intersection.Vertices';
end


% function checkVertexAngles(vertex_angles)
%   angle_diff = pwintz.math.angleDistance(vertex_angles)
% end