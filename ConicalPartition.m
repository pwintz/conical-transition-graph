classdef ConicalPartition < handle
  %%% ConicalPartition - Create a class for defining a partition of R^n into cones. 
  %%% To keep things simple, we are starting with 2D and 3D cases. 

  
  % ` runtests TestConicalPartition

  properties(Constant, GetAccess = {?matlab.unittest.TestCase}) % Define class constants.
    % Define the smallest angle between vertex angles before two vertices are merged. 
    ANGLE_TOL (1, 1) double = 1e-5;

    % Define the minimum distance between vertices for them to be unique.
    VERTEX_TOL (1, 1) double = 1e-3;
  end % End of class constants block

  properties(SetAccess = immutable) % Define public properties.
    % The dimension of the system considered
    dimension (1,1); 
    
    % ⋘──────── Counts ────────⋙
    
    % Number rays defining the facets of cones.
    n_rays         (1, 1); 
    % Number of (n-1) dimensional boundaries between cones (doesn't count boundaries between cones that only touch at a corner).
    n_boundaries   (1, 1);
    % Number of cones in the partition.
    n_cones        (1, 1);

    % Index arrays
    vertex_indices     (1, :); % 1:n_vertices
    cone_indices       (1, :); % 1:n_cones 
    boundary_indices   (1, :); % 1:n_boundaries
    ray_indices        (1, :); % 1:n_rays
    
    
    % The vertices array contains one vertex per column. Each vertex except for the origin is a unit vector. All vertices are unique.
    vertices       (:, :) double; 

    % ⋘──────── Vertices and Cones ────────⋙

    % Get the vertices without the origin.
    rays           (:, :) double; 
    
    % The cones cell array contains cone per entry.
    % cones          ConvexPolyhedralCone;
    cones          (1, :); % Array of ConvexPolyhedralCone.

    % Mesh.
    triangulation;

    % For convenience, we define the origin vector.
    origin         (:, 1) double; 
  end

  
  properties(SetAccess = immutable, GetAccess={?matlab.unittest.TestCase}) % Private properties

    
    % Number vertices in mesh (one for each ray, plus one for the origin.)
    n_vertices     (1, 1); 


    % Subsets of vertex_indices
    origin_index       (1, 1);
    ray_indices_in_vertices (1, :); % vertex_indices without the origin index. 
    
     
    
    mesh_graph; % An undirected mesh_graph of the vertices in the ConicalPartion with edges between if they both border the same cone.
    % cone_graph;   % An undirected mesh_graph of the cones in the ConicalPartion with edges between if they share a vertex (beside zero).
    vertex_table;
    cone_table;
    boundary_table;
  end

  methods(Static)

    function conical_partition = fromNumSlices(n_slices)
      % Partition R^2 usings a "watermelon slicing" method using equal angles for all cones.
      arguments(Input)
        n_slices (1,1) double;
      end % End of Input arguments block
      
      arguments(Output)
        conical_partition ConicalPartition;
      end % End of Output arguments block
    
      vertex_angles = pwintz.arrays.range(0, 2*pi, "n_values", n_slices, "includeEnd", false);
      conical_partition = ConicalPartition.fromAngles2D(vertex_angles);

      ConicalPartition.validateAngles(vertex_angles);
    end % End of function

    function conical_partition = fromAngles2D(vertex_angles)
      % Partition R^2 usings a "watermelon slicing" method using equal angles for all cones.
      arguments(Input)
        vertex_angles (1, :) double;
      end % End of Input arguments block
      
      % ⋘────────── Sort Angles and Remove Duplicates ──────────⋙
      vertex_angles = mod(vertex_angles, 2*pi); % Map to [0, 2*pi)
      vertex_angles = uniquetol(vertex_angles, ConicalPartition.ANGLE_TOL);
      % this.vertex_angles = sort(vertex_angles);
      
      % ⋘──────── Generate vertices ────────⋙
      origin = [0; 0];
      vertices = [origin, pwintz.math.angle2UnitVector(vertex_angles)];

      conical_partition = ConicalPartition(vertices);
    end % End of function

    function conical_partition = fromUVSphere3D(options)
      % Create a grid of points that covers the unit sphere.
      arguments(Input)
        options.nLinesOfLongitude = 10;
        options.nLinesOfLatitude  = 10;
      end % End of Input arguments block
      
      vertices = sphereGrid3d(nLinesOfLongitude=options.nLinesOfLongitude, nLinesOfLatitude=options.nLinesOfLatitude);
      conical_partition = ConicalPartition(vertices);
    end

  end % End static methods block

  methods

    % Constructor
    function this = ConicalPartition(vertices)
      arguments(Input)
        vertices (:, :) double;
      end % End of Input arguments block
      this.dimension = size(vertices, 1);

      this.origin = zeros(this.dimension, 1);
      is_zero_column = vecnorm(vertices) < 1e-12;
      rays = vertices(:, ~is_zero_column);
      rays = pwintz.arrays.normalizeColumns(rays);
      rays = pwintz.arrays.uniqueColumns(rays, tolerance=ConicalPartition.VERTEX_TOL);
      vertices = [rays, this.origin];

      % ╭────────────────────────────────────────────────╮
      % │             Generate Triangulation             │
      % ╰────────────────────────────────────────────────╯
      triang = constructTriangulation(vertices);

      % Update vertices to match the order in triang;
      vertices = triang.Points';
      this.vertices = vertices;
      
      n_vertices = size(triang.Points, 1);
      n_rays     = n_vertices - 1;

      switch this.dimension
        case 2
          % Get the indices of vertices in the boundary.  
          [~,bnd_vert_ndxs] = boundaryshape(triang);
          
          % Construct edges as sequential boundary vertices in the polygon produced by boundaryshape.
          bnd_edges = [bnd_vert_ndxs, circshift(bnd_vert_ndxs, 1)];
          
          % Each line in ConnectivityList is a cone.
          n_cones      = size(triang.ConnectivityList, 1); 
          % Get the number of (n-1)-dimensional boundaries between cones (i.e., does not include boundaries where only corners touch).
          n_boundaries = size(bnd_edges, 1); 
        case 3
          % Creating a boundary triangulation only works in 3D, where the boudnary is 2D, since MATLAB does not allow a 1D triangulation.
          [external_facets_vertex_ndxs, external_points] = triang.freeBoundary();
          boundary_triangulation = triangulation(external_facets_vertex_ndxs, external_points);
          
          % Each line in ConnectivityList is a cone.
          n_cones      = size(boundary_triangulation.ConnectivityList, 1);
          % (n-1)-dimensional boundaries between cones (i.e., does not include boundaries where only corners touch).
          n_boundaries = size(boundary_triangulation.edges(), 1); 
          
          % ⋘──────── Construct boudnary graph edges ────────⋙
          % Get a list of all edges on the boudnary of the mesh. The indices are relative to the "boundary_triangulation" mesh's indexing.
          bnd_edges = boundary_triangulation.edges();
          
          % Convert from the indices of "boundary_triangulation" to "triag"
          bnd_vert_ndxs = triang.nearestNeighbor(boundary_triangulation.Points);
          bnd_edges = bnd_vert_ndxs(bnd_edges);
        otherwise
          error("Unexpected case: %s.", this.dimension);
      end
      
      % Repackage the "triangulation" object into our own structure to make it compatible with existing code.
      mesh_graph = struct( ...
        "n_nodes", n_vertices, ...
        "n_edges", n_boundaries, ...
        "n_faces", n_cones, ...
        "nodes", triang.Points, ...
        "edges", triang.edges(), ...
        "faces", triang.ConnectivityList ...
      );

      rays        = pwintz.arrays.mapRows(@(pt_ndxs) {triang.Points(pt_ndxs, :)'}, triang.ConnectivityList);
      hspace_reps = pwintz.arrays.mapRows(@(ray_cell)  HalfspaceRepresentation.fromConicalHull(ray_cell{1}), rays);
      % hspace_reps = cell(n_cones, 1);
      % for cone_ndx = 1:n_cones
      %   vertex_ndxs = triang.ConnectivityList(cone_ndx, :);
      %   vertices    = triang.Points(vertex_ndxs, :)';
      %   hspace_reps{cone_ndx, 1} = HalfspaceRepresentation.fromConicalHull(vertices);
      % end
      % hspace_reps = [hspace_reps{:}]';

      cone_centers = triang.incenter();
      pwintz.strings.format("triang.Points %z\nthis.triangulation.ConnectivityList %z\nrays %z\nhspace_reps: %z %s\ncone_centers %z %s", triang.Points, triang.ConnectivityList, rays, hspace_reps, class(hspace_reps), cone_centers, class(cone_centers))
      
      % ⋘──────── Generate a mesh_graph that stores the relationships between vertices ────────⋙

      % ⋘──────── Save the origin's index and all of the other indices ────────⋙
      this.origin_index             =  pwintz.arrays.findRowIn(this.origin', triang.Points);
      ray_indices_in_vertices = setdiff(1:n_vertices, this.origin_index); 
      
      rays = triang.Points(ray_indices_in_vertices, :)';
      this.ray_indices = 1:n_rays; 
      this.rays = rays;
      this.ray_indices_in_vertices = ray_indices_in_vertices;

      % Check that the origin is excluded from the rays.
      pwintz.assertions.assertNotColumnIn(this.origin, rays, columnName="origin");

      % ╭────────────────────────────────────────────────╮
      % │             Construct Vertex Table             │
      % ╰────────────────────────────────────────────────╯

      pwintz.strings.format("Constructing vertex table for %d vertices.", n_vertices);
      cones_adjacent_to_node = triang.vertexAttachments();

      % Convert the indices to column vectors.
      cones_adjacent_to_node = cellfun(@(x) transpose(x), cones_adjacent_to_node, "UniformOutput", false);

%       triangulation_edges = triang.edges();
%       % edges_adjacent_to_node = 
%       edges_adjacent_to_node = cell(n_vertices, 1);
%       nodes_adjacent_to_node = cell(n_vertices, 1);
%       
% %       for i_vert = 1:n_vertices
%         % faces_adjacent_to_node{i_vert} = find(sum(mesh_graph.faces' == i_vert))';
% 
%         % Find all the rows in "triangulation_edges" that contains "i_vert". 
%         [vert_edge_ndxs, ~] = find(i_vert == triangulation_edges);
% 
%         % Find the endpoints of all edges connected to this 
%         verts_in_connected_edges = unique(triangulation_edges(vert_edge_ndxs, :));
% 
%         % Get the edges, from the edge indices.
%         edges_adjacent_to_node{i_vert} = triangulation_edges(vert_edge_ndxs, :);
% 
%         % Remove the current vertex.
%         nodes_adjacent_to_node{i_vert} = setdiff(verts_in_connected_edges, i_vert);
%       end

      vertex_table = struct2table(struct( ...
          "index", num2cell(1:n_vertices)', ... % Index is redundant, but makes reading the table visually easier.
          "position", num2cell(triang.Points', 1)', ...
          "adjacent_cones_ndxs", cones_adjacent_to_node ...
          ...."connected_mesh_edges_ndxs", edges_adjacent_to_node, ... % TODO Remove edges that are on the "outside" of the triangulation. This was not being used anywhere.
          ..."adjacent_node_ndxs", nodes_adjacent_to_node  ... % TODO: Update to "vertex"
      ));


      % ╭──────────────────────────────────────────────────╮
      % │             Construct Boundary Table             │
      % ╰──────────────────────────────────────────────────╯
      % Create a table that stores each boundary between cones that has dimension (n-1). 
      % Each entry should have
      % * The edge index from the mesh.
      % * The boundary cone object.
      % * The indices of the cones that are touching it non-trivially (should be only two).
      % * The indices of the vertices?

      pwintz.strings.format("Constructing boundary table for %d boundaries.", n_boundaries);
      cone_neighbors = triang.neighbors();
      % The "triangulation/neighbors()" function produces a 'nan' for each facet that has an external face (not touching any other triangles or tetrahedrons). We expect every facet to have an external face, so they should all have a nan entry. MATLAB puts the nan entries in the first column.
      assert(all( isnan(cone_neighbors(:, 1))), "We expect the first column to only contain NaN values.");
      cone_neighbors = cone_neighbors(:, 2:end);
      assert(all(~isnan(cone_neighbors), 'all'), "We expect to remove all of the nan entries by dropping the first column.");

      % ⋘──────── Construct an array that gives the original vertex index of each vertex in the boundary triangulation ────────⋙
      % original_ndxs_of_vertex_ndxs = triang.nearestNeighbor(boundary_triangulation.Points);
      
      % Get a list of all edges on the boudnary of the mesh. The indices are relative to the "boundary_triangulation" mesh's indexing.
%       bnd_edges = boundary_triangulation.edges();
% 
%       % Convert from the indices of "boundary_triangulation" to "triag"
%       bnd_edges = original_ndxs_of_vertex_ndxs(bnd_edges);

      % ⋘──────── Construct ConvexPolyhdralCones for each boundary ────────⋙
      bnd_cones = cell(n_boundaries, 1);
      for i_bnd = 1:n_boundaries
        edge = bnd_edges(i_bnd, :);
        bnd_cone_verts = triang.Points(edge, :)';
        bnd_cones{i_bnd} = ConvexPolyhedralCone.fromRays(bnd_cone_verts);
      end

      % Get the indices of the cones adjacent to each boundary edge.
      switch this.dimension
        case 2
          indices_of_adjacent_cones = triang.vertexAttachments(bnd_vert_ndxs);
        case 3
          indices_of_adjacent_cones = triang.edgeAttachments(bnd_edges);
        otherwise
          error("Unexpected case: %s.", this.dimension);
      end
      
      boundary_table = struct2table(struct(...
        "index", num2cell(1:n_boundaries)', ...
        "edge", num2cell(bnd_edges, 2),  ...
        "indices_of_adjacent_cones", indices_of_adjacent_cones,  ...
        "boundary_cone", bnd_cones ...
      ));

      % ╭──────────────────────────────────────────────╮
      % │             Construct Cone Table             │
      % ╰──────────────────────────────────────────────╯
      % Create a table that stores each cone in the conical partition. Each entry should have
      % * The indices of the vertices from the mesh that define the rays.
      % * The cone object.
      % * Adjacent cones.

      cone_centers = triang.incenter();

      pwintz.strings.format("Constructing cone table for %d cones.", n_cones);
      cones = cell(n_cones, 1);
      for i_cone = 1:n_cones
        % Get the rays defining each cone (and the origin)
        vertex_ndxs    = triang.ConnectivityList(i_cone, :);
        ray_ndxs       = setdiff(vertex_ndxs, this.origin_index);
        rays           = triang.Points(ray_ndxs, :)';
        pwintz.assertions.assertNumRows(rays, this.dimension);
        assert(~isempty(rays));
        cones{i_cone} = ConvexPolyhedralCone.fromRays(rays);
        assert(cones{i_cone}.n_rays == size(rays, 2));
        assert(cones{i_cone}.n_rays > 0);
      end

      cone_table = struct2table(struct( ...
          "cone_index", num2cell(1:n_cones)', ... % Index is redundant, but makes reading the table visually easier.
          "cone", cones, ... % Index is redundant, but makes reading the table visually easier.
          "vertex_indices", num2cell(triang.ConnectivityList, 2), ...
          "center_vector", num2cell(cone_centers, 2), ...
          "neighbor_cones_ndxs", num2cell(cone_neighbors, 2) ...
          ...."connected_mesh_edges_ndxs", edges_adjacent_to_node, ... % TODO Remove edges that are on the "outside" of the triangulation. This was not being used anywhere.
          ..."adjacent_node_ndxs", nodes_adjacent_to_node  ... % TODO: Update to "vertex"
      ));


      % Convert the cell array of cones to an array.
      this.cones = [cones{:}];

      
      this.mesh_graph     = mesh_graph;
      this.vertex_table   = vertex_table;
      this.boundary_table = boundary_table;
      this.cone_table     = cone_table;

      this.triangulation  = triang;

      % Counts.
      this.n_vertices     = n_vertices;
      this.n_rays         = n_rays;
      this.n_cones        = n_cones;
      this.n_boundaries   = n_boundaries;

      % Indices.
      this.vertex_indices   = 1:n_vertices;
      this.cone_indices     = 1:n_cones;
      this.boundary_indices = 1:n_boundaries;

      % ╭──────────────────────────────────────╮
      % │             Check values             │
      % ╰──────────────────────────────────────╯
      % We want the vertices and points in the triangulation to have the same order.
      pwintz.assertions.assertEqual(this.triangulation.Points', this.vertices);

      pwintz.assertions.assertSize(this.rays, [this.dimension, this.n_rays]);
      pwintz.assertions.assertSize(this.vertices, [this.dimension, this.n_vertices]);
      % pwintz.assertions.assertSize(this., [this.dimension, this.n_rays]);

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
      pwintz.strings.format("Finished constructing ConicalPartition: %s", this)
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
    
    function v = getVertex(this, vert_index)
      v = this.vertices(:, vert_index);
    end % End of function

    function ray = getRay(this, ray_index)
      assert(all(ray_index <= this.n_rays), "ray_index=%s must not be more than n_rays=%d", mat2str(ray_index), this.n_rays)
      ray = this.rays(:, ray_index);
    end % End of function
    
    function cone = getCone(this, index)
      % The cones property is an array of ConvexPolyhedralCones. This function returns a portion of that array with the same shape as index.
      
      cone = this.cones(index);
      pwintz.assertions.assertSize(cone, size(index));
    end % End of function

    function [bnd_cone, adjacent_cone_ndxs, ray_indices] = getBoundary(this, bnd_ndx)
      % The cones property is an array of ConvexPolyhedralCones. This function returns a portion of that array with the same shape as bnd_ndx.

      bnd_cone           = this.boundary_table.boundary_cone(bnd_ndx, :);
      adjacent_cone_ndxs = this.boundary_table.indices_of_adjacent_cones(bnd_ndx, :);
      ray_indices        = this.boundary_table.edges(bnd_ndx, :);

      % bnd = reshape(bnd, size(bnd_ndx));
      pwintz.assertions.assertSize(bnd_cone, size(bnd_ndx));
    end % End of function
    
    % ? unneeded? function [v, ndx] = getVertexFromAngle(this, angle)
    % ? unneeded?   % Find a list of vertices that match the given list of angles.
    % ? unneeded?   ndx  = find(pwintz.math.angleDistance(angle, this.vertex_angles) < ConicalPartition.ANGLE_TOL);
    % ? unneeded?   % ndx = this.normalizeVertexIndex(ndx);
    % ? unneeded?   v = this.vertices(ndx);
    % ? unneeded?   error("We need to get the index from among the vertices not counting the origin.");
    % ? unneeded? end % End of function

    % ╭─────────────────────────────────────────────────────╮
    % │             Vertex mesh_graph Structure             │
    % ╰─────────────────────────────────────────────────────╯
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

      adjacent_vertex_ndxs = grAdjacentNodes(this.mesh_graph.edges, v_ndx)';
      assert(all(adjacent_vertex_ndxs <= this.n_vertices), "adjacent_vertex_ndxs must be in {1, 2, ..., n_vertices}");
      assert(isrow(adjacent_vertex_ndxs), 'Output of getVerticesAdjacentToVertex must be a row vector');
    end % End of function


    function adjacent_cone_ndxs = getConesAdjacentToRay(this, ray_ndx)
      arguments(Input)
        this; 
        ray_ndx (1, 1) double {pwintz.validators.mustBeIndexScalar};
      end % End of Input arguments block
      arguments(Output)
        adjacent_cone_ndxs (1, :) double {pwintz.validators.mustBeIndexVector};
      end % End of Output arguments block

      assert(ismember(ray_ndx, this.ray_indices), "The index %d is not a ray index (meaning it is th origin).")

      vertex_ndx = this.ray_indices_in_vertices(ray_ndx);
      adjacent_cone_ndxs = this.getConesAdjacentToVertex(vertex_ndx);

      assert(numel(adjacent_cone_ndxs) ~= this.n_cones, "Too many cones were returned (%d). It seems that this ray '%s' is actually the origin.", numel(adjacent_cone_ndxs), mat2str(this.getVertex(vertex_ndx)));
    end

    function adjacent_cone_ndxs = getConesAdjacentToVertex(this, v_ndx)
      arguments(Input)
        this; 
        v_ndx (1, 1) double {pwintz.validators.mustBeIndexScalar};
      end % End of Input arguments block
      arguments(Output)
        adjacent_cone_ndxs (1, :) double {pwintz.validators.mustBeIndexVector};
      end % End of Output arguments block
      
      % !! Don't use this.triangulation.vertexAttachments to get the adjacent cones, since it only works in 3D.
      % x adjacent_cone_ndxs = this.triangulation.vertexAttachments(v_ndx);

      % this.triangulation.Points
      % this.boundary_table.indices_of_adjacent_cones
      % this.vertex_table.adjacent_cones_ndxs
      adjacent_cone_ndxs = this.vertex_table.adjacent_cones_ndxs(v_ndx, :);
      adjacent_cone_ndxs = cell2mat(adjacent_cone_ndxs);

      
      % adjacent_cone_ndxs = this.vertex_table.adjacent_cones_ndxs{v_ndx, :}';

      % ! We don't support multiple v_ndxs anymore. 
      % % When multiple vertex indices are given in v_ndx, its possible that they 
      % % share some adjacent cones, so we use 'unique' to filter duplicates. 
      % adjacent_cone_ndxs = unique(adjacent_cone_ndxs);
      
      % Sanity check.
      assert(all(adjacent_cone_ndxs <= this.n_cones), "adjacent_cone_ndxs must be in {1, 2, ..., n_cones}");
    end % End of function

    function result = areVerticesAdjacent(this, v1_ndx, v2_ndx)
      result = ismember(v2_ndx, grAdjacentNodes(this.mesh_graph.edges, v1_ndx));
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
      
      adjacent_vertex_ndxs = this.mesh_graph.faces(cone_ndx, :);

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
%       this.vertex_table.adjacent_cones_ndxs
%       this.vertex_table.adjacent_cones_ndxs(v_ndx, :)
% 
%       neighborConesNdxOffsets = [-1, 0]; % For 2D case
% 
%       adjacent_cones_ndxs = this.normalizeConeIndex(v_ndx + neighborConesNdxOffsets);
%     end


    function [intersecting_cone_ndxs] = getConesIntersectingArc(this, start_angle, end_angle)
      % Given start and end angles, this function construct an arc from the start angle to the end angle in the counter-clockwise direction. The angle swept by this arc must be less than or equal to pi. 
      % This function returns all of the cone indices that the arc intersects non-trivially. 
      % The angles are shifted slightly toward the center of the arc to introduce a small numerical margin.
      
      assert(this.dimension == 2, "This function only works in 2D.");
      assert(isscalar(start_angle), ...
        "start_angle must have shape %s but instead had shape %s", ...
        mat2str(size([1, 1])), mat2str(size(start_angle)));
      assert(isscalar(end_angle), ...
        "end_angle must have shape %s but instead had shape %s", ...
        mat2str(size([1, 1])), mat2str(size(end_angle)));
        
      start_angle = mod(start_angle, 2*pi);
      end_angle   = mod(end_angle, 2*pi);
      % assert(0 <= start_angle, "Must have 0 <= start_angle");
      % assert(0 <= start_angle, "Must have start_angle < 2*pi");
      % assert(0 <= end_angle, "Must have 0 <= end_angle");
      % assert(end_angle <= 2*pi, "Must have end_angle <= 2*pi");

      % disp("");
      % pwintz.strings.format("=== getConesIntersectingArc(start_angle=%.2f, end_angle=%.2f) ===", start_angle, end_angle)
      ray_angles = pwintz.math.atan2(this.rays);
      ray_angles = mod(ray_angles, 2*pi);
      ray_angles = sort(ray_angles);
      
      cone_midpoint_angles = ray_angles + pwintz.math.angleDiffCCW(ray_angles)/2;

      % Find the last ray that has an angle <= the start angle
      pwintz.strings.format("ray_angles: %.4g", ray_angles)
      last_ray_not_after_start_angle = find(ray_angles <= start_angle, 1, "last");
      first_ray_not_before_end_angle = find(ray_angles >= end_angle, 1, "first");
      if isempty(first_ray_not_before_end_angle)
        % If the angles of all of the rays is less than the end_angle, such then we loop around to use the first ray.
        first_ray_not_before_end_angle = 1;
      end

      % Find new start and end angles that are rounded to the nearest ray angle. 
      new_start_angle = ray_angles(last_ray_not_after_start_angle);
      new_end_angle   = ray_angles(first_ray_not_before_end_angle);

      % ?DEBUG pwintz.strings.format("    start_angle=%.3g,     end_angle=%.3g",     start_angle,     end_angle)
      % ?DEBUG pwintz.strings.format("new_start_angle=%.3g, new_end_angle=%.3g", new_start_angle, new_end_angle)

      % Since the angles loop around, we will either have a pattern of which rays are in the angles:
      %     [0, 0, 0, 1, 1, 1, 0]  % if start_ray_ndx < end_ray_ndx
      % or 
      %     [1, 1, 1, 0, 0, 1, 1] % if start_ray_ndx > end_ray_ndx
      if start_angle == end_angle 
        is_cone_midpoint_between_start_and_end = true(size(ray_angles));
      elseif new_start_angle < new_end_angle
        % ?DEBUG pwintz.strings.format("Case: start_angle = %.3g < end_angle = %.3g", start_angle, end_angle);
        is_cone_midpoint_between_start_and_end = new_start_angle <= cone_midpoint_angles & cone_midpoint_angles <= new_end_angle;
      else
        % ?DEBUG pwintz.strings.format("Case: start_angle = %.3g > end_angle = %.3g", start_angle, end_angle);
        is_cone_midpoint_between_start_and_end = new_start_angle <= cone_midpoint_angles | cone_midpoint_angles <= new_end_angle;
      end

      assert(any(is_cone_midpoint_between_start_and_end), "None of the angles where found to be in the range.");

      % midpoint_angle_indices
      
      midpoint_angles_in_arc = cone_midpoint_angles(is_cone_midpoint_between_start_and_end);

      % Construct vectors in the direction of the midpoints, scaled down to be within the polytope.
      midpoint_vectors = 0.1*pwintz.math.angle2UnitVector(midpoint_angles_in_arc);
      intersecting_cone_ndxs = this.triangulation.pointLocation(midpoint_vectors');
      return      

%       % Since "end_angle" can be numerically less than "start_angle" (e.g., if start_angle = 3 and end_angle = 0), we compute the angle difference between them in the counter-clockwise direction.  
%       arc_angle = pwintz.math.angleDiffCCW([start_angle, end_angle], index=1);
%       assert(arc_angle <= pi, "The angle from start_angle=%.2g to end_angle=%.2g in the counter-clockwise direction was %.2g but must be <= pi", start_angle, end_angle, arc_angle);
% 
%       % We nudge the angles slightly to avoid including cones that only touch on the boundary.
%       start_theta = start_angle             + ConicalPartition.ANGLE_TOL / 2;
%       end_theta   = start_angle + arc_angle - ConicalPartition.ANGLE_TOL / 2;
%       theta = linspace(start_theta, end_theta, 3);
%       points = [cos(theta); sin(theta)];
% 
%       % Generate the arc.
%       arc = Polytope.fromConvexHull(points);
%       intersecting_cone_ndxs = this.getConeIntersections(arc);
    end

    function cone_ndx_between_angles = getConesBetweenVerticesFromAngles(this, start_angle, end_angle)
      assert(this.dimension == 2, "This function only works in 2D.");
      assert(all(size(start_angle) == [1, 1]), ...
        "start_angle must have shape %s but instead had shape %s", mat2str(size([1, 1])), mat2str(size(start_angle)));
      assert(all(size(end_angle) == [1, 1]), ...
        "end_angle must have shape %s but instead had shape %s", mat2str(size([1, 1])), mat2str(size(end_angle)));

      % Bump the starting point slightly to the side so that we ensure it is only in one cone.
      start_point    = pwintz.math.angle2UnitVector(start_angle + ConicalPartition.ANGLE_TOL/2);
      end_point      = pwintz.math.angle2UnitVector(start_angle - ConicalPartition.ANGLE_TOL/2);
      start_cone_ndx = this.getConesContainingPoint(start_point);
      end_cone_ndx   = this.getConesContainingPoint(end_point);

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
      %   % for face = this.mesh_graph.faces
      %   %   if 
      %   % end
      % 
      %   adjacent_cones_ndxs = this.vertex_table.adjacent_cones_ndxs;
      %   is_eond_ndx_between_angles = zeros(size(this.mesh_graph.faces));
      %   for v_ndx = start_vertex_ndx:(end_vertex_ndx - 1)
      %     cone_ndx_in_set = intersect(adjacent_cones_ndxs{v_ndx}, adjacent_cones_ndxs{v_ndx+1});
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


    % function result = getInflowBoundaries(this, cone_ndx)
    %   arguments(Input)
    %     this;
    %     v (:, 1) double;
    %     bnd_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
    %     cone_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
    %   end % End of Input arguments block
    %   arguments(Output)
    %     result (1, 1) logical;
    %   end % End of Output arguments block
    %   cone = this.getCone(cone_ndx);
    %   bnd  = this.getBoundary(bnd_ndx);
    % end

%     function bnd_cone = doesVectorPointFromBoundaryIntoCone(this, v, bnd_ndx, cone_ndx)
%       arguments(Input)
%         this;
%         v (:, 1) double;
%         bnd_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
%         cone_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
%       end % End of Input arguments block
%       arguments(Output)
%         result (1, 1) logical;
%       end % End of Output arguments block
%       cone = this.getCone(cone_ndx);
%       bnd  = this.getBoundary(bnd_ndx);
% 
%       
%       
%     end % End of function

    function [normal_vectors, bnd_vertex_ndxs] = getConeNormals(this, cone_ndx) 
      % Find the vectors a_1 and a_2 such that the given cone equals {x | <x, a1> \geq 0, <x, a2> \geq 0}.
      % this.getCone(cone_ndx);

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

    function sphere_over_approx_poly = getUnitSphereOverApproximationInCone(this, cone_ndx)
      sphere_over_approx_poly = this.getCone(cone_ndx).getUnitSphereOverApproximation();
    end

    function middle_vector = getConeMiddleVector(this, cone_ndx) % Convenience function. Useful for plotting.
      % [~, adjacent_vertices] = this.getVerticesAdjacentToCone(cone_ndx);
      if nargin() == 1
        cone_ndx = this.cone_indices;
      end
      % this.cone_table.center_vector
      middle_vector = this.cone_table.center_vector(cone_ndx, :)';
      % middle_vector = mean(adjacent_vertices, 2);

      middle_vector = pwintz.arrays.normalizeColumns(middle_vector);
    end

    % ╭──────────────────────────────────────────╮
    % │             Set Computations             │
    % ╰──────────────────────────────────────────╯

    function cone_ndxs = getConesContainingPoint(this, x)
      % Find all of the cones that contain x. 
      % assert(this.dimension == 2, 'Only implemented for 2 dimensions');
      assert(iscolumn(x) && numel(x) == 2);

      % If x is at the origin, then it is in all of the cones.
      if  norm(x) == 0
        cone_ndxs = this.cone_indices;
        return
      end % End of if " norm(x) == 0" block

      % Rescale x so that it within the polytope defined b 
      x = 0.1 * x / norm(x);
      cone_ndxs = this.triangulation.pointLocation(x');
      return
      % cone_ndx = tri

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
      % Find the 2 closest vertices to x.
      k = 2;
      [~, closest_vertex_ndxs] = mink(vecnorm(this.vertices - x), k);
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
        intersection = cone.intersection(polyhedron);
        intersections{cone_ndx} = intersection;
        does_cone_intersect(cone_ndx) = ~isempty(intersection);
        fprintf(".");
      end
      fprintf("\n");

      % Convert from a logical array to indices.
      intersecting_cone_ndxs = find(does_cone_intersect);
    end % End of function

    % ╭───────────────────────────────────────────────────╮
    % │             Overload Built-in methods             │
    % ╰───────────────────────────────────────────────────╯
    function disp(this)
      builtin("disp", this);
      if all(size(this) == [1, 1])
        disp("Vertex Graph");
        disp(this.mesh_graph);
        disp("Vertex Positions");
        disp(this.mesh_graph.nodes);
        disp("Edges between vertices");
        disp(this.mesh_graph.edges);
        disp("Cone Vertices");
        disp(this.mesh_graph.faces);
        disp("Vertex table (Private)");
        disp(this.vertex_table);
      end
    end % End of function

    function string_rep = char(this)
      % This method defines the string representation that is 
      % inserted when using "%s" in a string format, such as 
      % fprintf("%s", polyhedron) or sprintf("%s", polyhedron); 
      % It does not change the "disp" output. 
      arguments(Output)
        string_rep char; % Ensure output is cast to char, even if you create a string.
      end
      string_rep = sprintf("%s(dimension=%d, n_vertices=%d, n_rays=%d, n_boundaries=%d, n_cones=%d)", class(this), this.dimension, this.n_vertices, this.n_rays, this.n_boundaries, this.n_cones);
    end
    % function disp(this)
    %   fprintf('%s\n', this);
    % end % End of function

    function plot(this, varargin)
      colors = {0.5*[1 1 1]}; % Gray.
      switch this.dimension
        case 2
          for cone_ndx = this.cone_indices
            cone = this.getCone(cone_ndx);
            color = colors{mod(cone_ndx, numel(colors)) + 1}; 
            this.plotCone(cone_ndx, "FaceColor", color, "FaceAlpha", 0.2);
          end
        case 3
          tetramesh(this.triangulation, 'FaceAlpha', 0.2);
        otherwise
          error("Unexpected case: dimension=%s.", this.dimension);
      end
      
      % ⋘──────── plot middle vectors ────────⋙
      % cone_middle_vector = this.conical_partition.getConeMiddleVector(reachable_cone_ndx);
      % pwintz.plots.plotVector2(cone_middle_vector, plotArgs={':k', "ShowArrowHead", false, "HandleVisibility", "off"});

      % ⋘──────── Plot the normal vectors ────────⋙
      % [normal_vectors, bnd_vertex_ndxs] = this.conical_partition.getConeNormals(reachable_cone_ndx);
      % normal_at_i_vertex = normal_vectors(:, bnd_vertex_ndxs == i_vertex_ndx);
      % pwintz.plots.plotVector2(0.5*i_vertex, 0.1*normal_at_i_vertex, plotArgs={"Color", [0.4, 0.0, 0.7], "HandleVisibility", "off"});
    end % End of function

    function plotCone(this, cone_ndx, varargin)
      switch this.dimension
        case 2
          T = this.triangulation.ConnectivityList(cone_ndx, :);
          x = this.triangulation.Points(:, 1);
          y = this.triangulation.Points(:, 2);
          z = 0*x;
          % The "triplot" function draws outlines of triangulation in 2D, but doesn't seem to allow filling the curves, and also has different options than "tetramesh", which makes it annoying to handle both, so we don't use it.
          % x triplot(T, x, y, "displayName", sprintf("Cone %d", cone_ndx));

          % Instead, "trimesh" draws filled triangles for triangulation, but only works for 3D meshes, so we add third component z=0.
          trimesh(T, x, y, z, "EdgeColor", "black", "displayName", sprintf("Cone %d", cone_ndx), varargin{:});
        case 3
          tetramesh(this.triangulation.ConnectivityList(cone_ndx, :), this.triangulation.Points, 'FaceAlpha', 0.2, "displayName", sprintf("Cone %d", cone_ndx), varargin{:});
        otherwise
          error("Unexpected case: %s.", this.dimension);
      end
      
    end % End of function

    function plotBoundary(this, bnd_ndx)
      bnd_cone = this.boundary_table.boundary_cone(bnd_ndx, :);
      adjacent_cone_ndxs = this.boundary_table.indices_of_adjacent_cones(bnd_ndx, :);
      
      label = sprintf("Cone %d \\cap Cone %d", adjacent_cone_ndxs(1), adjacent_cone_ndxs(2));

      verts = [bnd_cone.rays, this.origin];

      switch this.dimension
        case 2
          plot(verts(1,:), verts(2,:), 'r');
        case 3
          color = rand(1, 3);
          fill3(verts(1,:), verts(2,:), verts(3,:), color, "DisplayName", label);
        otherwise
          error("Unexpected case: %s.", this.dimension);
      end
    end % End of function

    function plotOverapproximationOfUnitSphereInCone(this, cone_ndx)
      cone = this.getCone(cone_ndx);
      sphere_approx = cone.getUnitSphereOverApproximation();
      % label = sprintf("Cone %d \\cap Cone %d", adjacent_cone_ndxs(1), adjacent_cone_ndxs(2));
      sphere_approx.plot();
    %   adjacent_cone_ndxs = this.boundary_table.indices_of_adjacent_cones(bnd_ndx, :);
    % 
    %   verts = [bnd_cone.rays, this.origin];
    % 
    %   switch this.dimension
    %     case 2
    %       plot(verts(1,:), verts(2,:), 'r');
    %     case 3
    %       color = rand(1, 3);
    %       fill3(verts(1,:), verts(2,:), verts(3,:), color, "DisplayName", label);
    %     otherwise
    %       error("Unexpected case: %s.", this.dimension);
    %   end
    end % End of function

    function plotOverapproximationOfUnitSphereInBoundary(this, bnd_ndx)
      bnd_cone = this.boundary_table.boundary_cone(bnd_ndx, :);
      sphere_approx = bnd_cone.getUnitSphereOverApproximation();
      % label = sprintf("Cone %d \\cap Cone %d", adjacent_cone_ndxs(1), adjacent_cone_ndxs(2));
      sphere_approx.plot();
    %   adjacent_cone_ndxs = this.boundary_table.indices_of_adjacent_cones(bnd_ndx, :);
    % 
    %   verts = [bnd_cone.rays, this.origin];
    % 
    %   switch this.dimension
    %     case 2
    %       plot(verts(1,:), verts(2,:), 'r');
    %     case 3
    %       color = rand(1, 3);
    %       fill3(verts(1,:), verts(2,:), verts(3,:), color, "DisplayName", label);
    %     otherwise
    %       error("Unexpected case: %s.", this.dimension);
    %   end
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


function tri = constructTriangulation(vertices)
  dim = size(vertices, 1);
  origin = zeros(dim, 1);
  del_tri = delaunayTriangulation(vertices');
  
  % ⋘──────── Cleanup triangulation ────────⋙
  % The initial triangulation may have some triangles that don't include the origin, so we find the boundary of the triangulation and for each triangle in the boudnary, we add the origin manually.
  [vertex_ndxs, rays] = del_tri.freeBoundary();

  n_cones = size(vertex_ndxs, 1);

  % ⋘──────── Add origin to vertices ────────⋙
  % Insert the origin at the beginning
  points = [origin'; rays];
  % Put the origin's index (1) into the vertex_ndxs array. Increment all of the existing indices by +1 because we each point down to fit the origin in the first index.
  vertex_ndxs = [ones(n_cones, 1), vertex_ndxs+1];

  % % Check that all of the vertices are used.
  % pwintz.assertions.assertAllAreMembers(1:size(points), unique(vertex_ndxs));
  tri = triangulation(vertex_ndxs, points);

end