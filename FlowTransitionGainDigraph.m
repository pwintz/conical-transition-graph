classdef FlowTransitionGainDigraph < TransitionGainDigraph
  % ` runtests TestFlowTransitionGainDigraph

  properties % Define instance variables.
    flow_set_cone_ndxs (1, :) {pwintz.validators.mustBeIndexVector};
    flow_map_matrix    (:, :) {pwintz.validators.mustBeSquare};
    can_flow_from_ray_into_cone (:, :) logical;
    can_flow_from_cone_to_ray   (:, :) logical;
    % The (i,j) entry of directly_reachable_sets_from_ray_through_cone contains 
    % the set that is reachable from vertex i by flowing through cone j.
    directly_reachable_sets_from_ray_through_cone   (:,:) cell; % ConvexPolyhedron; 
    restricted_reachable_sets_from_unit_sphere_in_cone (1,:) cell; % ConvexPolyhedron;
  end % End of properties block

  methods
    function this = FlowTransitionGainDigraph(conical_partition, flow_set_cone_ndxs, flow_map_matrix)
      this@TransitionGainDigraph(conical_partition);
      
      this.flow_set_cone_ndxs = flow_set_cone_ndxs;
      this.flow_map_matrix = flow_map_matrix;
      % cone_indices = conical_partition.cone_indices;
      % is_cone_in_flow_set = ismember(cone_indices, flow_set_cone_ndxs);

      % ╭───────────────────────────────────────────────────────────────╮
      % │             Cone-Vertex Reachability Calculations             │
      % ╰───────────────────────────────────────────────────────────────╯
      % ⋘──────── Construct a list of cones that can be reached from each vertex ────────⋙
      % ⋘──────── Construct a list of vertices that can be reached from each cone ────────⋙
      can_flow_from_ray_into_cone = zeros(conical_partition.n_rays, conical_partition.n_cones);
      can_flow_from_cone_to_ray   = zeros(conical_partition.n_cones, conical_partition.n_rays);
      
      for i_ray_ndx = conical_partition.ray_indices
        ray = conical_partition.getRay(i_ray_ndx);
        flow_dir = flow_map_matrix * ray; % Flow direction at "ray".

        % If the cone is not in the flow set, then it's impossible to flow into it.
        adjacent_flow_set_cone_ndxs = intersect(conical_partition.getConesAdjacentToRay(i_ray_ndx), flow_set_cone_ndxs);
        for adjacent_cone_ndx = adjacent_flow_set_cone_ndxs
          adj_cone = conical_partition.getCone(adjacent_cone_ndx);
          
          if adj_cone.atXDoesVPointInward(ray, flow_dir)
            can_flow_from_ray_into_cone(i_ray_ndx, adjacent_cone_ndx) = true;
          else
            can_flow_from_cone_to_ray(i_ray_ndx, adjacent_cone_ndx) = true;
          end
      
      % x     % Compute <Ac*v, n> where \dot x = Ac*x is the system's flow, v is a vertex in the boundary of the cone, and n is the normal vector to the cone.
      % x     [normal_vectors, bnd_vertex_ndxs] = conical_partition.getConeNormals(adjacent_cone_ndx);
      % x     normal_at_ray = normal_vectors(:, bnd_vertex_ndxs == i_ray_ndx);
      % x     % Compute \dot x = Ax for x = ray.
      % x     flow_at_ray = flow_map_matrix * ray; 
      % x     pwintz.assertions.assertSize(normal_at_ray, [2, 1]);
      % x 
      % x     % Sanity check.
      % x     assert(abs(dot(normal_at_ray, ray)) < 1e-6, 'The normal and vector must be orthogonal.');
      
          % % If <Ac*v, n> >= 0, then it is possible to flow into the cone from v.
          % if dot(normal_at_ray, flow_at_ray) >= 0
          %   can_flow_from_ray_into_cone(i_ray_ndx, adjacent_cone_ndx) = true;
          %   % We use a (min/max) gain of 1.0 when flowing from a vertex to a cone because the distance flowed to reach the cone is zero.  
          %   this.addEdgeFromVertexToCone(i_ray_ndx, adjacent_cone_ndx, 1.0, 1.0);
          % else
          %   can_flow_from_cone_to_ray(adjacent_cone_ndx, i_ray_ndx) = true;
          %   this.addEdgeFromConeToVertex(adjacent_cone_ndx, i_ray_ndx, 0, 1);
          % end
        end % End of "adjacent_cone_ndx" for loop.
      end % End "i_ray_ndx" for
      
      % ╭────────────────────────────────────────────────────────────────╮
      % │             Direction(s) of flows at each boundary             │
      % ╰────────────────────────────────────────────────────────────────╯
      % ! In 3D, we need to determine flow across2D  boundaries instead of vertices, but in 2D the boundaries are vertices.
%       can_flow_from_bnd_into_cone = zeros(conical_partition.n_boundaries, conical_partition.n_cones);
%       can_flow_from_cone_into_bnd = zeros(conical_partition.n_cones, conical_partition.n_boundaries);
%       for i_bnd_ndx = conical_partition.boundary_indices
%         [bnd, adjacent_cone_ndxs, ray_indices] = conical_partition.getBoundary(i_bnd_ndx);
% 
%         can_flow_from_ray_into_cone(i_ray_ndx, adjacent_cone_ndx) = true;
%         can_flow_from_cone_to_ray(i_ray_ndx, adjacent_cone_ndx) = true;
% 
%         can_flow_from_bnd_into_cone(i_bnd_ndx)
%         can_flow_from_cone_into_bnd(i_bnd_ndx)
% 
%         rays = bnd.rays;
%         for ray = rays 
% 
%         end
%       end
      
      % can_flow_from_ray_into_cone
      % can_flow_from_ray_into_cone(conical_partition.origin_index, :)
      % assert(all(~can_flow_from_ray_into_cone(conical_partition.origin_index, :)), "No cones can be reached from the origin.");
      this.can_flow_from_ray_into_cone = can_flow_from_ray_into_cone;
      this.can_flow_from_cone_to_ray   = can_flow_from_cone_to_ray;
      
      % ╭──────────────────────────────────────────────────────────────────╮
      % │             Intra-cone Reachability from Unit Sphere             │
      % ╰──────────────────────────────────────────────────────────────────╯
      % Compute the region in each cone C_i that can be reached by a solution to \dot x = A_c * C_i that starts in unit sphere is always in C_i.
      restricted_reachable_sets_from_unit_sphere_in_cone = cell(conical_partition.n_cones, 1);
      for cone_ndx = this.flow_set_cone_ndxs
        sphere_overapproximation = conical_partition.getUnitSphereOverApproximationInCone(cone_ndx);
        [restricted_reachable_set, ~] = this.computeReachableSet(cone_ndx, sphere_overapproximation);
        assert(~isempty(restricted_reachable_set))
        restricted_reachable_sets_from_unit_sphere_in_cone{cone_ndx} = restricted_reachable_set;
      end
      this.restricted_reachable_sets_from_unit_sphere_in_cone = restricted_reachable_sets_from_unit_sphere_in_cone;
      
      % ╭────────────────────────────────────────────────────────────────────────────────────────────╮
      % │             Contstruct Flow Graph for Intra-cone Reachability between Vertices             │
      % ╰────────────────────────────────────────────────────────────────────────────────────────────╯
      directly_reachable_sets_from_ray_through_cone = cell(conical_partition.n_rays, conical_partition.n_cones);
      % directly_reachable_sets_from_ray_through_cone = ConvexPolyhedron.arrayOfEmpty(conical_partition.n_vertices, conical_partition.n_cones);
      
      for ray_ndx = conical_partition.ray_indices
        cone_ndx = find(this.can_flow_from_ray_into_cone(ray_ndx, :));
        if isempty(cone_ndx)
          continue
        end
        % assert(isscalar(cone_ndx), "Is it possible that a flow only points in from one vertex? Guess we'll see! cone_ndx=%s", mat2str(cone_ndx))

        ray    = conical_partition.getRay(ray_ndx);
        
        reach_set = this.computeReachableSet(cone_ndx, Polytope.fromPoint(ray));
        if isempty(reach_set)
          error('The reach set from from ray_ndx=%d in cone_ndx=%d is empty: %s', ray_ndx, cone_ndx, reach_set);
        end
        directly_reachable_sets_from_ray_through_cone{ray_ndx, cone_ndx} = reach_set;
      

        assert(reach_set.isBounded(), "Reachable set is unbounded.");
        reach_set_vertices = reach_set.vertices;
        initial_ray_ndx_in_reach_set_vertices = pwintz.arrays.findColumnIn(ray, reach_set_vertices, tolerance=1e-3, verbose=false);

        % ⋘──────── Get the reachable set vertices excluding the initial ray ────────⋙
        
        reach_set_vertices = reach_set_vertices(:, setdiff(1:end, initial_ray_ndx_in_reach_set_vertices));
      
        ray_norms = vecnorm(reach_set_vertices); % Used to compute gains.
        normalized_reach_set_vertices = nrv(reach_set_vertices);
        for bnd_vertex_ndx = conical_partition.getVerticesAdjacentToCone(cone_ndx)
          bnd_vertex = conical_partition.getVertex(bnd_vertex_ndx);
          
          col_ndxs = pwintz.arrays.findColumnIn(bnd_vertex, normalized_reach_set_vertices, tolerance=1e-5, verbose=false);

          if ~isempty(col_ndxs)
            gains = ray_norms(col_ndxs);
            max_gain = max(gains);
            min_gain = min(gains);
            % Add transition from the non-origin "vertex_ndx" to the adjacent vertex "bnd_vertex_ndx" that can be reached by a single interval of flow.
            this.addEdgeFromVertexToVertex(ray_ndx, bnd_vertex_ndx, min_gain, max_gain);
          end
        end
        % this.flow_graph = flow_graph;
      
      end % End for loop
      this.directly_reachable_sets_from_ray_through_cone = directly_reachable_sets_from_ray_through_cone;

      assert(~isempty(directly_reachable_sets_from_ray_through_cone));

    end % End of constructor

    % ╭────────────────────────────────────────╮
    % │             Reachable Sets             │
    % ╰────────────────────────────────────────╯
    function [reachable_set_in_cone, reachable_set] = computeReachableSet(this, cone_ndx, P0)
      % Compute the cone that is reachable from P0 when flowing according to \dot x \in A*C_i, where C_i is given cone .
      % The "derivative cone" is represented using vertices on the unit sphere, so we scale it up so that we get a set the sufficiently represents "P0 + A*C_i", at least up to a large radius.
      C_i  = this.conical_partition.getCone(cone_ndx);
      AC_i = this.flow_map_matrix * C_i;
      
      % !! This code is not robust. We are using convex polytopes to represent the cones and just scaling them to make them hopefully "sufficiently close to infinite", but when we make the number large, the intersection calculations break down for numerical reasons.
      
      reachable_set         = P0 + AC_i;
      reachable_set_in_cone = intersection(reachable_set, C_i);
    end

    function transition_graph = contractEdgesToBetweenCones(this, start_cone_ndxs, end_cone_ndxs)
      arguments(Input)
        this FlowTransitionGainDigraph;
        start_cone_ndxs (1, :) double {mustBeInteger,mustBePositive};
        end_cone_ndxs   (1, :) double {mustBeInteger,mustBePositive};
      end % End of Input arguments block
      arguments(Output)
        transition_graph (1, 1) TransitionGainDigraph;
      end % End of Output arguments block

      % ` Psuedocode:
      % `
      % ` new_graph = new TransitionGainGraph. 
      % ` for each start index and each end index
      % `    search for a path from the start cone to the end cone in the form 
      % `            [cone] -> [boundary] -> ... -> [boundary] -> [cone].
      % `    for each path found, 
      % `       compute the max and min gains
      % `       add edge to new_graph
      % `    end
      % ` end
      % Check arguments.
           
      transition_graph = TransitionGainDigraph(this.conical_partition);

      start_cone_ndxs = intersect(start_cone_ndxs, this.flow_set_cone_ndxs);
      end_cone_ndxs   = intersect(end_cone_ndxs,   this.flow_set_cone_ndxs);

      % Check that the cone indices are OK.
      pwintz.assertions.assertAllAreMembers(start_cone_ndxs, this.conical_partition.cone_indices);
      pwintz.assertions.assertAllAreMembers(end_cone_ndxs,   this.conical_partition.cone_indices);
      
      % Construct the subgraph from that only includes edges between vertices.
      edges_from_verts_to_verts = this.getEdgesFromVerticesToVertices();
      vertex_subgraph = digraph(edges_from_verts_to_verts, this.gains_digraph.Nodes);


      edges_not_from_verts_to_verts = [this.getEdgesFromConesToVertices(); this.getEdgesFromVerticesToCones()];
      faces_subgraph = digraph(edges_not_from_verts_to_verts, this.gains_digraph.Nodes);
      
      % % We will also use edges from vertices to cones...
      % edges_from_verts_to_cones = this.getEdgesFromVerticesToCones();
      % % ...and from cones to verts.
      % edges_from_cones_to_verts = this.getEdgesFromConesToVertices();

      % Get all of the vertices that are next to the start and end cones.
      vertices_adjacent_to_start_cones = this.conical_partition.getVerticesAdjacentToCone(start_cone_ndxs);
      vertices_adjacent_to_end_cones   = this.conical_partition.getVerticesAdjacentToCone(end_cone_ndxs);

      for start_vertex_ndx = vertices_adjacent_to_start_cones
        start_node_ndx = this.node_ndxs_of_vertices(start_vertex_ndx);
        for end_vertex_ndx = vertices_adjacent_to_end_cones
          end_node_ndx = this.node_ndxs_of_vertices(end_vertex_ndx);
          
          % fprintf('Finding all paths from vertex %3d to vertex %3d.\n', start_vertex_ndx, end_vertex_ndx);
          % !! Previously we skipped cases where the vertices are equal, but this actually an important case, since we want to include paths like 
          % `   [cone] -> [vertex] -> [cone]
          % !! which only has path lengths of length zero. 
          [paths_nodes, paths_edges] = vertex_subgraph.allpaths(start_vertex_ndx, end_vertex_ndx);
          if isempty(paths_nodes)
            % No paths found.
            fprintf('Did not find any paths from vertex %3d (node %d) to vertex %3d (node %d).\n', start_vertex_ndx, start_node_ndx, end_vertex_ndx, end_node_ndx);
            pwintz.assertions.assertNotEqual(start_vertex_ndx, end_vertex_ndx);
            continue
          end

          for i_path = pwintz.arrays.cellVectorIndices(paths_nodes)
          % for i_path = 1:numel(paths_nodes)
            path_nodes = paths_nodes{i_path};
            path_edges = paths_edges{i_path};
            
            assert(path_nodes(1)   == start_node_ndx, 'path_nodes=%s must start at start_node_ndx=%d', mat2str(path_nodes), start_node_ndx);
            assert(path_nodes(end) == end_node_ndx, 'path_nodes=%s must end at end_node_ndx=%d', mat2str(path_nodes), end_node_ndx);
            % fprintf('Found path for start_cone_ndx=%d, end_cone_ndx=%d, start_vertex_ndx=%d, end_vertex_ndx=%d, i_path=%d\n', start_cone_ndx, end_cone_ndx, start_vertex_ndx, end_vertex_ndx, i_path);
            fprintf('Found path from vertex %3d (node %d) to vertex %3d (node %d): %s.\n', start_vertex_ndx, start_node_ndx, end_vertex_ndx, end_node_ndx, mat2str(path_nodes));

            % Now that we have a path through the vertices, we need to add the start and end edges that connect to the adjacent cones.
            for start_cone_ndx = start_cone_ndxs
              start_cone_node_ndx = this.node_ndxs_of_cones(start_cone_ndx);
              start_edge = this.getEdgesFromConeToVertex(start_cone_ndx, start_vertex_ndx);
              if isempty(start_edge) 
                fprintf('\tNo start_edge from start_cone_ndx=%d (node %d) to start_vertex_ndx=%d. Continuing to next start cone.\n', start_cone_ndx, start_cone_node_ndx, start_vertex_ndx);
                continue
              else
                fprintf('\tFound a start edge from start_cone_ndx=%d (node %d) to start_vertex_ndx=%d. Continuing to next start cone.\n', start_cone_ndx, start_cone_node_ndx, start_vertex_ndx);
              end
              for end_cone_ndx = end_cone_ndxs
                end_cone_node_ndx   = this.node_ndxs_of_cones(end_cone_ndx);
                end_edge = this.getEdgesFromVertexToCone(end_vertex_ndx, end_cone_ndx);
                if isempty(end_edge) 
                  fprintf('\t\tNo end_edge from end_vertex_ndx=%d to end_cone_ndx=%d (node %d). Continuing to next end cone.\n', end_vertex_ndx, end_cone_ndx, end_cone_node_ndx);
                  continue
                else
                  fprintf('\t\tFound a complete path! edges from end_vertex_ndx=%d to end_cone_ndx=%d (node %d). Continuing to next end cone.\n', end_vertex_ndx, end_cone_ndx, end_cone_node_ndx);

                  % Check: There should only be at most one edge between a cone and a vertex.
                  pwintz.assertions.assertNumRows(start_edge, 1);
                  pwintz.assertions.assertNumRows(end_edge,   1);

                  % ⋘──────── Construct the complete path from cone to cone ────────⋙
                  path = [
                    start_edge;
                    this.gains_digraph.Edges(path_edges, :); % <- path_edges may be empty.
                    end_edge
                  ];
                  
                  % Compute the cumulative gains as the product of all of the edge gains.
                  min_gain = prod(path.MinGain);
                  max_gain = prod(path.MaxGain);

                  % ⋘──────── Add the path between cone  ────────⋙
                  transition_graph.addEdgeFromConeToCone(start_cone_ndx, end_cone_ndx, min_gain, max_gain);
                end
                
                % start_node_ndx 
                % start_cone_ndxs
                % edges_from_verts_to_cones.EndNodes == [start_node_ndx, start_cone_ndx] 
                % end_node_ndx
                % edges_from_cones_to_verts
              end
            end
            % start_cone_node_ndxs
            % end_cone_node_ndxs
            % edges_to_search_for = [start_vertex_ndx*ones(size(start_cone_ndxs')), start_cone_node_ndxs']
            % ismember(edges_from_verts_to_cones.EndNodes, edges_to_search_for, 'rows')


          end % End of "path" for loop.

        end % End of for loop
      end % End of for loop

    end % End of function



%   end % End of methods block.
% 
%   methods(Access = {?matlab.unittest.TestCase}) % Private block, accessible by unit tests.
  end % End private methods block

end % End of class