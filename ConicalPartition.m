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

  properties(SetAccess = immutable) % Define public properties.
    % The dimension of the system considered
    dimension (1,1); 
    
    % Counts
    n_vertices     (1, 1); % The numbers of state vertices and derivative vertices are equal.
    n_cones        (1, 1); % The numbers of state vertices and derivative vertices are equal.
   
    % Index arrays
    vertex_indices (1, :);
    cone_indices   (1, :);
    origin_index   (1, 1);
    nonorigin_vertex_indices (1, :);

    % ⋘──────── Vertices and Cones ────────⋙
    % The vertices array contains one vertex per column. Each vertex except for the origin is a unit vector. All vertices are unique.
    vertices       (2, :) double; 
    % The cones cell array contains cone per entry. Each vertex except for the origin is a unit vector. All vertices are unique.
    cones          (1, :) cell;

    % For convenience, we define the origin vector.
    origin         (2, 1) double; 
  end

  % TODO: Make private after debugging:, GetAccess=private)
  properties(SetAccess = immutable, Hidden) % Private properties
    vertex_multi_index_size (1,1);
    vertex_angles      (1, :) double;     
    
    vertex_graph; % An undirected vertex_graph of the vertices in the ConicalPartion with edges between if they both border the same cone.
    % cone_graph;   % An undirected vertex_graph of the cones in the ConicalPartion with edges between if they share a vertex (beside zero).
    vertex_table;
    % cone_table;
  end

  methods(Static)
    function test(varargin)% Define convenience function for running tests.
      TestConicalPartition.runTests(varargin{:});
    end % End of function

    function conical_partition = fromNumSlices(n_slices)
      arguments(Input)
        n_slices (1,1) double;
      end % End of Input arguments block
      
      arguments(Output)
        conical_partition ConicalPartition;
      end % End of Output arguments block
    
      vertex_angles = pwintz.arrays.range(0, 2*pi, "n_values", n_slices, "includeEnd", false);
      conical_partition = ConicalPartition(vertex_angles);
    end % End of function

    % ! Move to ConicalAbstraction? function conical_partition = fromAngles(A, derivative_vertex_angles, extra_state_vertex_angles)
    % ! Move to ConicalAbstraction?   % p = inputParser;
    % ! Move to ConicalAbstraction?   % addRequired(p, "A", @(A) isnumeric(A) && ismatrix(A) && issquare(A)); % Positional argument.
    % ! Move to ConicalAbstraction?   % addParameter(p, "nSlices", 8, @(x) isreal(x) && rem(x,1)==0 && x > 1); % Name-value parameter.;
    % ! Move to ConicalAbstraction?   % parse(p, varargin{:});
    % ! Move to ConicalAbstraction?
    % ! Move to ConicalAbstraction?   if isempty(extra_state_vertex_angles)
    % ! Move to ConicalAbstraction?     combined_derivative_vertex_angles = derivative_vertex_angles;
    % ! Move to ConicalAbstraction?   else
    % ! Move to ConicalAbstraction?     extra_state_vertices = pwintz.math.angle2UnitVector(extra_state_vertex_angles);
    % ! Move to ConicalAbstraction?     extra_derivative_vertex_angles = pwintz.math.atan2(A * extra_state_vertices);
    % ! Move to ConicalAbstraction?
    % ! Move to ConicalAbstraction?     % Combine the given derivative vertex angles with the ones coming from the extra state vertex angles. We don't need to sort them here because that is done in the ConicalPartition constructor.
    % ! Move to ConicalAbstraction?     combined_derivative_vertex_angles = [derivative_vertex_angles, extra_derivative_vertex_angles];
    % ! Move to ConicalAbstraction?   end
    % ! Move to ConicalAbstraction?
    % ! Move to ConicalAbstraction?   % vertices_table = struct2table(struct(...
    % ! Move to ConicalAbstraction?   %     "state_vertex",      state_vertices, ...
    % ! Move to ConicalAbstraction?   %     "derivative_vertex", derivative_vertices, ...
    % ! Move to ConicalAbstraction?   %     "state_angle",       state_vertex_angles, ...
    % ! Move to ConicalAbstraction?   %     "derivative_angle",  derivative_vertex_angles  ...
    % ! Move to ConicalAbstraction?   % ))
    % ! Move to ConicalAbstraction?   conical_partition = ConicalPartition(A, combined_derivative_vertex_angles);
    % ! Move to ConicalAbstraction? end % End of function
  end % End static methods block

  methods

    % Constructor
    function this = ConicalPartition(vertex_angles)
      arguments(Input)
        vertex_angles (1, :) double;
      end % End of Input arguments block
      
      this.dimension = 2;
      this.vertex_multi_index_size = 1; % For 2D Case. 

      % ⋘────────── Sort Angles and Remove Duplicates ──────────⋙
      vertex_angles = mod(vertex_angles, 2*pi); % Map to [0, 2*pi)
      vertex_angles = uniquetol(vertex_angles, ConicalPartition.ANGLE_TOL);
      % this.vertex_angles = sort(vertex_angles);

      % ⋘──────── Generate vertices ────────⋙
      origin = [0; 0];
      this.vertices = [origin, pwintz.math.angle2UnitVector(vertex_angles)];
      this.origin = origin;

      % ⋘──────── Generate a vertex_graph that stores the relationships between vertices ────────⋙
      [nodes, edges] = gabrielGraph(this.vertices');
      this.origin_index =  pwintz.arrays.findRowIn([0, 0], nodes);
      not_origin_ndxs = find(1:size(nodes, 1) ~= this.origin_index);
      this.nonorigin_vertex_indices = not_origin_ndxs;

      % ⋘──────── Construct faces ────────⋙
      faces_ndx_col_1 = not_origin_ndxs';
      faces_ndx_col_2 = circshift(not_origin_ndxs,1)';
      faces_ndx_col_3 = this.origin_index * ones(size(not_origin_ndxs'));
      faces = [faces_ndx_col_1, faces_ndx_col_2, faces_ndx_col_3];

      % node_to_adjacent_face_ndxs = dictionary(k1,v1,...,kN,vN);

      vertex_graph = struct( ...
        "n_nodes", size(nodes, 1), ...
        "n_edges", size(edges, 1), ...
        "n_faces", size(faces, 1), ...
        "nodes", nodes, ...
        "edges", edges, ...
        "faces", faces ...
      );

      % ╭────────────────────────────────────────────────╮
      % │             Construct Vertex Table             │
      % ╰────────────────────────────────────────────────╯
      faces_adjacent_to_node = cell(vertex_graph.n_nodes, 1);
      edges_adjacent_to_node = cell(vertex_graph.n_nodes, 1);
      nodes_adjacent_to_node = cell(vertex_graph.n_nodes, 1);
      for i_node = 1:vertex_graph.n_nodes
        faces_adjacent_to_node{i_node} = find(sum(vertex_graph.faces' == i_node))';
        edges_adjacent_to_node{i_node} = grAdjacentEdges(vertex_graph.edges, i_node); 
        nodes_adjacent_to_node{i_node} = grAdjacentNodes(vertex_graph.edges, i_node);
        
        % [vertex_graph.edges(vertex_graph.edges(:, 1) == i_node, 2); 
                                          %  vertex_graph.edges(vertex_graph.edges(:, 2) == i_node, 1)];
      end
      vertex_table = struct2table(struct( ...
          "index", num2cell(1:vertex_graph.n_nodes)', ... % Index is redundant, but makes reading the table visually easier.
          "position", num2cell(nodes, 2), ...
          "adjacent_face_ndxs", faces_adjacent_to_node, ...
          "adjacent_edge_ndxs", edges_adjacent_to_node, ...
          "adjacent_node_ndxs", nodes_adjacent_to_node  ...
      ));

      this.vertex_graph = vertex_graph;
      this.vertex_table = vertex_table;


      this.n_vertices     = vertex_graph.n_nodes;
      this.n_cones        = vertex_graph.n_faces;
      this.vertex_indices = 1:vertex_graph.n_nodes;
      this.cone_indices   = 1:vertex_graph.n_faces;

      % ╭─────────────────────────────────────────────────╮
      % │             Construct list of Cones             │
      % ╰─────────────────────────────────────────────────╯
      cones = cell(this.n_cones, 1);
      for cone_ndx = this.cone_indices
        face_vertices = this.getVertex(this.vertex_graph.faces(cone_ndx, :));
        cones{cone_ndx} = ConvexPolyhedron.fromConvexHull(face_vertices);
      end
      this.cones = cones;

      % ╭──────────────────────────────────────────────╮
      % │             Construct Cone Table             │
      % ╰──────────────────────────────────────────────╯
      
      % this.cone_table = struct2table(struct(...
      %     "index", num2cell(1:vertex_graph.n_faces)', ...
      %     "vertex_indices", num2cell(faces, 2) ...
      % ));

% 
%       % ⋘──── Compute Derivative Vertex Angles of Extra State Vertex Angles ────⋙
%       for i = 1:numel(extra_state_vertex_angles)
%         extra_state_vertices(:,i)           = [cos(extra_state_vertex_angles(i)); sin(extra_state_vertex_angles(i))];
%         extra_derivative_vertex_angles(:,i) = pwintz.math.atan2(A * extra_state_vertices(:,i));
%       end
%
%       % ⋘──── Compute Property Values Derived from Given Values ────⋙
%       derivative_vertices = pwintz.math.angle2UnitVector(vertex_angles);
%       state_vertices      = A \ derivative_vertices;
%       state_vertex_angles = pwintz.math.atan2(state_vertices);
% 
%       % Remap from [-pi, pi] to [0, 2pi].
%       state_vertex_angles = mod(state_vertex_angles, 2*pi); 

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

      ConicalPartition.validateAngles(vertex_angles);
      % It is important that derivative_vertex_angles is sorted because we use that property in getDerivativeConeIndexContainingPoint.
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

    % ╭──────────────────────────────────╮
    % │             Indexing             │
    % ╰──────────────────────────────────╯
    % ? unneeded? function multi_indexArray = getVertexMultiIndices(this)
    % ? unneeded?   % vertexMultiIndices - Get an wide array where each column is the multi-index of a one vertex in this conical partition.
    % ? unneeded?   % The result can be used like 
    % ? unneeded?   % *
    % ? unneeded?   % *  for i = conical_partition.getVertexMultiIndices()
    % ? unneeded?   % *    v_i = conical_partition.getStateSpaceVertex(i)
    % ? unneeded?   % *  end
    % ? unneeded?   % *
    % ? unneeded?   assert(this.dimension == 2, "Only 2D implemented.");

    % ? unneeded?   % For the 2D case, each multi_indexing only has one entry.
    % ? unneeded?   % multi_indexArray = 1:this.n_vertices; %  <- For 2D.
    % ? unneeded?   % multi_indexArray = 1:this.n_slices;
    % ? unneeded?   multi_indexArray = this.vertex_indices;
    % ? unneeded? end % End of function
    
    function v = getVertex(this, index)
      v = this.vertices(:, index);
    end % End of function

    function cone = getConeVertices(this, index)
      cone = this.cones(index);
    end % End of function

    function cone = getCone(this, index)
      % The cones property is a cell array. This function returns a cell array with the same shape as index.
      cone = [this.cones{index}];
      cone = reshape(cone, size(index));
      pwintz.assertions.assertSize(cone, size(index));
    end % End of function
    
    % ? unneeded? function [v, ndx] = getVertexFromAngle(this, angle)
    % ? unneeded?   % Find a list of vertices that match the given list of angles.
    % ? unneeded?   ndx  = find(pwintz.math.angleDistance(angle, this.vertex_angles) < ConicalPartition.ANGLE_TOL);
    % ? unneeded?   % ndx = this.normalizeVertexIndex(ndx);
    % ? unneeded?   v = this.vertices(ndx);
    % ? unneeded?   error("We need to get the index from among the vertices not counting the origin.");
    % ? unneeded? end % End of function

    % ╭───────────────────────────────────────────────────────╮
    % │             Vertex vertex_graph Structure             │
    % ╰───────────────────────────────────────────────────────╯
    function adjacent_vertex_ndxs = getVerticesAdjacentToVertex(this, v_ndx)
      arguments(Input)
        this; 
        v_ndx (:, 1) double {pwintz.validators.mustBeIndexVector};
      end % End of Input arguments block
      arguments(Output)
        adjacent_vertex_ndxs (1, :) double {pwintz.validators.mustBeIndexVector};
      end % End of Output arguments block

      % !!! Using adjacent nodes works in 2D, but in higher dimensions, to find all of the vertices that are on the boundary of the same cone, we need to use a different approach.
      assert(this.dimension == 2);

      adjacent_vertex_ndxs = grAdjacentNodes(this.vertex_graph.edges, v_ndx)';
      assert(all(adjacent_vertex_ndxs <= this.n_vertices), "adjacent_vertex_ndxs must be in {1, 2, ..., n_vertices}");
      assert(isrow(adjacent_vertex_ndxs), 'Output of getVerticesAdjacentToVertex must be a row vector');
    end % End of function

    function adjacent_cone_ndxs = getConesAdjacentToVertex(this, v_ndx)
      arguments(Input)
        this; 
        v_ndx (:, 1) double {pwintz.validators.mustBeIndexVector};
      end % End of Input arguments block
      arguments(Output)
        adjacent_cone_ndxs (1, :) double {pwintz.validators.mustBeIndexVector};
      end % End of Output arguments block
      
      adjacent_cone_ndxs = this.vertex_table.adjacent_face_ndxs{v_ndx, :}';

      % When multiple vertex indices are given in v_ndx, its possible that they 
      % share some adjacent cones, so we use 'unique' to filter duplicates. 
      adjacent_cone_ndxs = unique(adjacent_cone_ndxs);
      
      % Sanity check.
      assert(all(adjacent_cone_ndxs <= this.n_cones), "adjacent_cone_ndxs must be in {1, 2, ..., n_cones}");
    end % End of function

    function result = areVerticesAdjacent(this, v1_ndx, v2_ndx)
      result = ismember(v2_ndx, grAdjacentNodes(this.vertex_graph.edges, v1_ndx));
    end % End of function

    function conjoining_cones_ndxs = getConesConjoiningTwoVertices(this, v1_ndx, v2_ndx)
      arguments(Input)
        this; 
        v1_ndx (1, 1) double {mustBeInteger,mustBePositive,mustBeNonempty,mustBeScalarOrEmpty};
        v2_ndx (1, 1) double {mustBeInteger,mustBePositive,mustBeNonempty,mustBeScalarOrEmpty};
      end % End of Input arguments block
      
      % getConesConjoiningTwoVertices - Get a list of cones (as indices) that are adjacent to two given vertices.
      % The result is wide array where each column is the multi-index of a cone.
      % v1_ndx
      % v2_ndx
      % assert(v1_ndx ~= v2_ndx, "The vertices passed to getConesConjoiningTwoVertices must be distinct. Instead we had v1_ndx=%d, v2_ndx=%d", v1_ndx, v2_ndx)
      cones_adjacent_to_v1 = this.getConesAdjacentToVertex(v1_ndx);
      cones_adjacent_to_v2 = this.getConesAdjacentToVertex(v2_ndx);
    
      % The conjoing cones are all cones that are adjacent to both vertices.
      conjoining_cones_ndxs = intersect(cones_adjacent_to_v1, cones_adjacent_to_v2);

      % if this.dimension == 2 && v1_ndx ~= this.origin_index && v2_ndx ~= this.origin_index && v1_ndx ~= v2_ndx
      %   assert(numel(conjoining_cones_ndxs) == 1);
      % end
    end

    function [adjacent_vertex_ndxs, adjacent_vertices] = getVerticesAdjacentToCone(this, cone_ndx, options)
      arguments(Input)
        this;
        cone_ndx (:, 1) {mustBeInteger,mustBePositive};
        options.includeOrigin logical = false;
      end % End of Input arguments block
      arguments(Output)
        adjacent_vertex_ndxs (1, :) {pwintz.validators.mustBeIndexVector};
        adjacent_vertices    (:, :) double;
      end % End of Output arguments block
      
      
      adjacent_vertex_ndxs = this.vertex_graph.faces(cone_ndx, :);

      % When multiple cone indices are given in cone_ndx, its possible that they 
      % share some adjacent vertices, so we use 'unique' to filter duplicates. 
      adjacent_vertex_ndxs = unique(adjacent_vertex_ndxs);

      % If includeOrigin=true, then delete the origin.
      if ~options.includeOrigin
        adjacent_vertex_ndxs = adjacent_vertex_ndxs(adjacent_vertex_ndxs ~= this.origin_index);
      end
      
      % If two output arguments, then find the vertices from their indices.
      if nargin() == 2
        adjacent_vertices = this.getVertex(adjacent_vertex_ndxs);
      end
    end % End of function

    % ╭──────────────────────────────────────────────────────────────────╮
    % │             Translating between vertices and cones             │
    % ╰──────────────────────────────────────────────────────────────────╯
    
%     function v_neighbors_ndxs = getVerticesAdjacentToVertex(this, multi_index)
%       % getVerticesAdjacentToVertex - Get an wide array containing where each column is the multi-index of a vertex that is adjacent to the given multi-index.
%       % The result can be used like 
%       % *                                                
%       % * for i = conical_partition.getVertexMultiIndices() 
%       % *   v_i = conical_partition.getStateSpaceVertex(i)         
%       % *   for j = conical_partition.getVerticesAdjacentToVertex(i) 
%       % *     v_j = conical_partition.getStateSpaceVertex(i)  
%       % *   end                                                   
%       % * end                                                   
%       % *                                                
% 
%       neighborKernel = [-1, +1]; % For 2D case
% 
%       v_neighbors_ndxs = this.normalizeVertexIndex(multi_index + neighborKernel);
%     end

%      function adjacent_cones_ndxs = getConesAdjacentToVertex(this, v_ndx)
%       % getConesAdjacentToVertex - Get a list of multi-indices for cones with the given vertex (defined by v_ndx) in their boundary.
%       % The result is wide array where each column is the multi-index of a cone.
%       v_ndx = this.normalizeVertexIndex(v_ndx);
%       
%       this.vertex_table.adjacent_face_ndxs
%       this.vertex_table.adjacent_face_ndxs(v_ndx, :)
% 
%       neighborConesNdxOffsets = [-1, 0]; % For 2D case
% 
%       adjacent_cones_ndxs = this.normalizeConeIndex(v_ndx + neighborConesNdxOffsets);
%     end


    function [intersecting_cone_ndxs, arc] = getConesIntersectingArc(this, start_angle, end_angle)
      % Given start and end angles, this function construct an arc from the start angle to the end angle in the counter-clockwise direction. The angle swept by this arc must be less than or equal to pi. 
      % This function returns all of the cone indices that the arc intersects non-trivially. 
      % The angles are shifted slightly toward the center of the arc to introduce a small numerical margin.
      
      assert(this.dimension == 2, "This function only works in 2D.");
      assert(isscalar(start_angle), ...
        "start_angle must have shape %s but instead had shape %s", mat2str(size([1, 1])), mat2str(size(start_angle)));
      assert(isscalar(end_angle), ...
        "end_angle must have shape %s but instead had shape %s", mat2str(size([1, 1])), mat2str(size(end_angle)));


      % Since "end_angle" can be numerically less than "start_angle" (e.g., if start_angle = 3 and end_angle = 0), we compute the angle difference between them in the counter-clockwise direction.  
      arc_angle = pwintz.math.angleDiffCCW([start_angle, end_angle], index=1);
      assert(arc_angle <= pi, "The angle from start_angle=%.2g to end_angle=%.2g in the counter-clockwise direction was %.2g but must be <= pi", start_angle, end_angle, arc_angle);

      % We nudge the angles slightly to avoid including cones that only touch on the boundary.
      start_theta = start_angle             + ConicalPartition.ANGLE_TOL / 2;
      end_theta   = start_angle + arc_angle - ConicalPartition.ANGLE_TOL / 2;
      theta = linspace(start_theta, end_theta, 3);
      points = [cos(theta); sin(theta)];

      % Generate the arc.
      arc = ConvexPolyhedron.fromConvexHull(points);
      intersecting_cone_ndxs = this.getConeIntersections(arc);

      % this.plotCones(intersecting_cone_ndxs, "")
    end

    
    function cone_ndx_between_angles = getConesBetweenVerticesFromAngles(this, start_angle, end_angle)
      assert(this.dimension == 2, "This function only works in 2D.");
      assert(all(size(start_angle) == [1, 1]), ...
        "start_angle must have shape %s but instead had shape %s", mat2str(size([1, 1])), mat2str(size(start_angle)));
      assert(all(size(end_angle) == [1, 1]), ...
        "end_angle must have shape %s but instead had shape %s", mat2str(size([1, 1])), mat2str(size(end_angle)));

      
      % Bump the starting point slightly to the side so that we ensure it is only in one cone.
      start_point = pwintz.math.angle2UnitVector(start_angle + ConicalPartition.ANGLE_TOL/2);
      end_point = pwintz.math.angle2UnitVector(start_angle - ConicalPartition.ANGLE_TOL/2);
      start_cone_ndx = this.getConesContainingPoint(start_point);
      end_cone_ndx = this.getConesContainingPoint(end_point);

      if end_cone_ndx >= start_cone_ndx
        cone_ndx_between_angles = start_cone_ndx : end_cone_ndx;
      else
        cone_ndx_between_angles = [1:end_cone_ndx, start_cone_ndx:this.n_cones];
      end

      %   start_cone_ndx:end_cone_ndx
      % 
      %   start_vertex_ndx = this.getVertexFromAngle(start_angle);
      %   end_vertex_ndx   = this.getVertexFromAngle(end_angle);
      %   % start_vertex_ndx  = find(pwintz.math.angleDistance(, ) < tol);
      %   % end_vertex_ndx    = find(pwintz.math.angleDistance(this.vertex_angles,   end_angle) < tol);
      %   assert(numel(start_vertex_ndx) > 0);
      %   assert(numel(  end_vertex_ndx) > 0);
      %   if numel(start_vertex_ndx) > 1
      %     error("Found multiple elements in this.vertex_angles=%s that match start_angle=%f, namely at indices=%s", mat2str(this.vertex_angles), start_angle, mat2str(start_vertex_ndx))
      %   end
      %   if numel(end_vertex_ndx) > 1
      %     error("Found multiple elements in this.vertex_angles=%s that match end_angle=%f, namely at indices=%s", mat2str(this.vertex_angles), end_angle, mat2str(end_vertex_ndx))
      %   end
      % 
      %   % for face = this.vertex_graph.faces
      %   %   if 
      %   % end
      % 
      %   adjacent_face_ndxs = this.vertex_table.adjacent_face_ndxs;
      %   is_eond_ndx_between_angles = zeros(size(this.vertex_graph.faces));
      %   for v_ndx = start_vertex_ndx:(end_vertex_ndx - 1)
      %     cone_ndx_in_set = intersect(adjacent_face_ndxs{v_ndx}, adjacent_face_ndxs{v_ndx+1});
      %     is_eond_ndx_between_angles(cone_ndx_in_set) = 1;
      %   end
      % 
      %   % assert(all(size(start_vertex_ndx) == [1, 1]), ...
      %   %   "start_vertex_ndx must be 1x1 but instead had shape %s", mat2str(size(start_vertex_ndx)));
      %   % assert(all(size(end_vertex_ndx) == [1, 1]), ...
      %   %   "end_vertex_ndx must be 1x1 but instead had shape %s", mat2str(size(end_vertex_ndx)));
      %   % assert(start_vertex_ndx ~= end_vertex_ndx);
      %   % 
      %   % if start_vertex_ndx < end_vertex_ndx
      %   %   cone_ndxs = start_vertex_ndx:(end_vertex_ndx - 1);
      %   % else
      %   %   % Handle the case where the indices loop past 0.
      %   %   cone_ndxs = [start_vertex_ndx:this.n_vertices, 1:(end_vertex_ndx - 1)];
      %   % end
    end

    function [normal_vectors, bnd_vertex_ndxs] = getConeNormals(this, cone_ndx) 
      % Find the vectors a_1 and a_2 such that the given cone equals {x | <x, a1> \geq 0, <x, a2> \geq 0}.
      [bnd_vertex_ndxs, bnd_vertices] = this.getVerticesAdjacentToCone(cone_ndx, includeOrigin=false);
      bnd_vertex_angles = pwintz.math.atan2(bnd_vertices);
      angle_counter_clockwise_diff = pwintz.math.angleDiffCCW(bnd_vertex_angles);
      is_cone_swept_ccw_from_v1_to_v2 = angle_counter_clockwise_diff(2) > angle_counter_clockwise_diff(1);
      if is_cone_swept_ccw_from_v1_to_v2
        % The normal vector of 
        warning('getConeNormals has not been tested for is_cone_swept_ccw_from_v1_to_v2=true.');
        normal_vector_angles = bnd_vertex_angles + [+pi/2, -pi/2];
      else
        normal_vector_angles = bnd_vertex_angles + [-pi/2, +pi/2];
      end
      normal_vectors = pwintz.math.angle2UnitVector(normal_vector_angles);
    end

    function max_angle = getConeMaxAngle(this, cone_ndx)
      % Compute the largest angle between boundary vertices for the given cone.
      % To find the largest angle, we use the fact that the angle theta between the vertices v_i and v_j (which are unit vectors) is given by cos(theta) = v_i'*v_j. The value of theta is maximized for the choice of (i,j) such that v_i'*v_j is minimized.
      assert(this.dimension == 2, 'Only tested for 2 dimensions, but it should work.');
      [~, adjacent_vertices] = this.getVerticesAdjacentToCone(cone_ndx);
      % Compute V'*V, where V = [v1, v2, ..., vp] is the matrix formed from all the cones boundary vertices. 
      dot_products_between_each_vertex = adjacent_vertices'*adjacent_vertices;
      min_dot_product = min(dot_products_between_each_vertex, [], 'all');
      max_angle = acos(min_dot_product);
    end

    function sphere_over_approx_poly = getUnitSphereOverApproximationInCone(this, cone_ndx)
      assert(this.dimension == 2, 'Only implemented for 2 dimensions');
      [~, adjacent_vertices] = this.getVerticesAdjacentToCone(cone_ndx);
      max_angle = this.getConeMaxAngle(cone_ndx);
      assert(max_angle < pi, 'max_angle must be less than pi.');
      % From geometry we have that convex hull of {v1, v2, r*v1, r*v2} contains the unit sphere if  .
      r = sec(max_angle / 2); 
      sphere_over_approx_poly = ConvexPolyhedron.fromConvexHull([adjacent_vertices, r*adjacent_vertices]);
    end

    function middle_vector = getConeMiddleVector(this, cone_ndx) % Convenience function. Useful for plotting.
      [~, adjacent_vertices] = this.getVerticesAdjacentToCone(cone_ndx);
      middle_vector = mean(adjacent_vertices, 2);
      middle_vector = middle_vector / norm(middle_vector);
    end

    % ╭──────────────────────────────────────────╮
    % │             Set Computations             │
    % ╰──────────────────────────────────────────╯

    function cone_ndxs = getConesContainingPoint(this, x)
      % Find all of the cones that contain x. 
      assert(this.dimension == 2, 'Only implemented for 2 dimensions');
      assert(iscolumn(x) && numel(x) == 2);

      % If x is at the origin, then it is in all of the cones.
      if  norm(x) == 0
        cone_ndxs = this.cone_indices;
        return
      end % End of if " norm(x) == 0" block

      % x_angle = pwintz.math.atan2(x);
      % x_angle = mod(x_angle, 2*pi); % Convert range to [0, 2pi] instead of [-pi, pi].
      % assert(x_angle >= 0);
      % 
      % % The values derivative_vertex_angles are sorted, by definition, so if we take the last index where the angle is greater than the vertex angle, then we'll have that it is the cone where x is between that vertex and the next vertex.
      % padded_vertex_angles = [this.vertex_angles, this.vertex_angles(1) + 2*pi];
      % 
      % % The index of the cone that x is in the index of the first angle that x_angle is larger than.
      % % !! Bug: This code fails to account for the fact that the indexing in the vertex angles is different than the vertices.
      % cone_ndxs = find(x_angle >= padded_vertex_angles, 1, 'last');
      % 
      % assert(~isempty(cone_ndxs), 'cone_ndxs must not be empty. ');

      % As an alternative method for finding the vertices, we can use the fact that the cone a point belongs to will have two vertices closer to the point than any other cone.
      x = x/norm(x);
      [~, closest_vertex_ndxs] = mink(vecnorm(this.vertices - x), 2);
      % cone_vertex_ndxs = [closest_vertex_ndxs, this.origin_index]; 
      alt_cone_ndxs =  this.getConesConjoiningTwoVertices(closest_vertex_ndxs(1), closest_vertex_ndxs(2));
      % assert(all(alt_cone_ndxs == cone_ndxs), "alt_cone_ndxs = %s, cone_ndxs = %s", mat2str(alt_cone_ndxs), mat2str(cone_ndxs))
      cone_ndxs = alt_cone_ndxs;
    end % End of function

    function [intersecting_cone_ndxs, does_cone_intersect, intersections] = getConeIntersections(this, polyhedron)
      % Get all of the cones that intersect with a given polyhedral.

      arguments(Input)
        this ConicalPartition;
        polyhedron ConvexPolyhedron; % <- We should make a new data type for polys.
      end
      arguments(Output)
        intersecting_cone_ndxs (1, :); 
        does_cone_intersect    (1, :) logical;
        intersections          (1, :) cell;
      end

      does_cone_intersect = zeros(size(this.cone_indices));
      intersections       = cell(size(this.cone_indices));

      for cone_ndx = this.cone_indices
        cone = this.getCone(cone_ndx);
        intersection  = 1e5*cone | polyhedron;
        intersections{cone_ndx} = intersection;
        does_cone_intersect(cone_ndx) = ~isempty(intersection);
      end

      % Convert from a logical array to indices.
      % TODO: when cone indices are multi-indices, we might need to convert the output of find, which gives linear indices, using ind2sub, or something like that.
      intersecting_cone_ndxs = find(does_cone_intersect);
    end % End of function

    % ╭───────────────────────────────────────────────────╮
    % │             Overload Built-in methods             │
    % ╰───────────────────────────────────────────────────╯
    function disp(this)
      builtin("disp", this);
      if all(size(this) == [1, 1])
        disp("Vertex Graph");
        disp(this.vertex_graph);
        disp("Vertex Positions");
        disp(this.vertex_graph.nodes);
        disp("Edges between vertices");
        disp(this.vertex_graph.edges);
        disp("Cone Vertices");
        disp(this.vertex_graph.faces);
        disp("Vertex table (Private)");
        disp(this.vertex_table);
      end
    end % End of function

    function plot(this, varargin)
      colors = {0.5*[1 1 1]}; % Gray.
      for cone_ndx = this.cone_indices
        color = colors{mod(cone_ndx, numel(colors)) + 1}; 
        plot(this.getCone(cone_ndx), "FaceColor", color, "FaceAlpha", 0.2);
      end

      % ⋘──────── plot middle vectors ────────⋙
      % cone_middle_vector = this.conical_partition.getConeMiddleVector(reachable_cone_ndx);
      % pwintz.plots.plotVector2(cone_middle_vector, plotArgs={':k', "ShowArrowHead", false, "HandleVisibility", "off"});

      % ⋘──────── Plot the normal vectors ────────⋙
      % [normal_vectors, bnd_vertex_ndxs] = this.conical_partition.getConeNormals(reachable_cone_ndx);
      % normal_at_i_vertex = normal_vectors(:, bnd_vertex_ndxs == i_vertex_ndx);
      % pwintz.plots.plotVector2(0.5*i_vertex, 0.1*normal_at_i_vertex, plotArgs={"Color", [0.4, 0.0, 0.7], "HandleVisibility", "off"});
    end % End of function

  end  % End of methods


  methods(Access = {?matlab.unittest.TestCase}) % Define private class methods.

%     function out = normalizeVertexIndex(this, multi_index)
%       % normalizeVertexIndex - Given a multi-index, normalize its entries by applying the modulo operation so that each entry is in {1, 2, ..., n_slices}.
%       assert(this.dimension == 2, "Only 2D implemented.");
% 
%       % Check that the multi-index matches the dimension of the system.
%       % assert(~isempty(multi_index), "multi_index must be nonempty.");
%       assert( ...
%         size(multi_index, 1) == this.vertex_multi_index_size, ...
%         "The size(multi_index, 1) = %d must be vertex_multi_index_size = %d.", size(multi_index, 1), this.vertex_multi_index_size ...
%       )
% 
%       % Compute the index of the vertex from the given argument, modulo the number of slices. Since MATLAB uses 1-based indexing, we need to subtract 1 before computing the modulo and then add it afterward.
%       out = mod(multi_index - 1, this.n_vertices) + 1;
%     end
% 
%     function out = normalizeConeIndex(this, index)
%       % normalizeVertexIndex - Given a multi-index, normalize its entries by applying the modulo operation so that each entry is in {1, 2, ..., n_slices}.
%       assert(this.dimension == 2, "Only 2D implemented.");
%   
%       % Check that the multi-index matches the dimension of the system.
%       % assert( ...
%       %   size(multi_index, 1) == this.cone_multi_index_size, ...
%       %   "The size(multi_index, 1) = %d must be vertex_multi_index_size = %d.", size(multi_index, 1), this.vertex_multi_index_size ...
%       % )
%   
%       % Compute the index of the vertex from the given argument, modulo the number of slices. Since MATLAB uses 1-based indexing, we need to subtract 1 before computing the modulo and then add it afterward.
%       out = mod(index - 1, this.n_cones) + 1;
%     end
  end% End of private methods

  methods(Static, Access = {?matlab.unittest.TestCase}) % Some static functions for internal use.
    
      function validateAngles(vertex_angles)
        % ⋘────────── Check range ──────────⋙
        if any(vertex_angles < 0 | vertex_angles >= 2*pi)
          out_of_range_ndxs = find(vertex_angles < 0 | vertex_angles >= 2*pi);
          error(...
            "ConicalPartition:AnglesOutOfRange", ...
            'The values of vertex_angles must in [0, 2*pi). The values %s are out of range at index(s) %s', ...
            mat2str(out_of_range_ndxs), ...
            mat2str(find(out_of_range_ndxs))...
          );
        end

        % ⋘────────── Check strictly increasing ──────────⋙
        if any(diff(vertex_angles) <= 0)
          error(...
            "ConicalPartition:AnglesNonincreasing", ...
            'The values of vertex_angles must be increasing. The angle decrease at index(s) %s', ...
            mat2str(find(diff(vertex_angles) > 0))...
          );
        end
        
        % ⋘────────── Check step size ──────────⋙
        % Check the step size of the angles in angle_steps as measured in the counter-clockwise direction. 
        angle_steps = pwintz.math.angleDiffCCW(vertex_angles);
        if any(angle_steps > pi)
          big_step_ndxs = find(angle_steps > pi);
          error(...
            "ConicalPartition:AngleStepTooBig", ...
            'The difference between consecutive angles in vertex_angles was %s > pi at index(s). ', ...
            mat2str(angle_steps(big_step_ndxs)), ...
            mat2str(big_step_ndxs) ...
          );
        end

        pwintz.assertions.assertUnique(vertex_angles, tolerance=ConicalPartition.ANGLE_TOL);

        % ╭───────────────────────────────────────────────────╮
        % │             Check vertex_angles             │
        % ╰───────────────────────────────────────────────────╯
        if nargin() == 2 % The second "state_angle_steps" can be omitted for testing purposes.
          % ⋘────────── Check length matches derivative_vertex_angles ──────────⋙
          % assert(all(size(derivative_vertex_angles) == size(state_vertex_angles)));

%           % ⋘────────── Check that there are no duplications ──────────⋙
%           [unique_angles, ndxs_of_unique_angles] = uniquetol(state_vertex_angles, ConicalPartition.ANGLE_TOL, 'OutputAllIndices', true);
%           if(numel(unique_angles) < numel(state_vertex_angles)) % If there were duplicates:
%             % Generate an informative error message...
%             is_entry_duplicate = cellfun(@(entry) numel(entry) > 1, ndxs_of_unique_angles');
%             duplicate_entries_str = cellfun(@(ndxs) sprintf("%f is duplicated at indices %s", state_vertex_angles(ndxs(1)), mat2str(ndxs)), ndxs_of_unique_angles(is_entry_duplicate))
% 
%             % Raise error.
%             error("There was a repeated entry in 'state_vertex_angles' (up to a small numerical tolerance). The list had %d items but only %d were unique. List of duplicates: \n\t%s.", numel(state_vertex_angles), numel(ndxs_of_unique_angles), join(duplicate_entries_str, "\n\t"));
%           end
%           
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
  
end % End of class