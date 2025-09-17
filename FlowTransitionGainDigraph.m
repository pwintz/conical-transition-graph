classdef FlowTransitionGainDigraph < TransitionGainDigraph
  % ` runtests TestFlowTransitionGainDigraph

  properties % Define instance variables.
    flow_set_cone_ndxs cell;
    flow_map_matrix    cell;
    % can_flow_from_vertex_into_cone (:, :) logical;
    % can_flow_from_cone_to_vertex   (:, :) logical;
    % The (i,j) entry of directly_reachable_sets_from_ray_through_cone contains 
    % the set that is reachable from vertex i by flowing through cone j.
    % directly_reachable_sets_from_ray_through_cone   (:,:) cell; % ConvexPolyhedron; 
    % restricted_reachable_sets_from_unit_sphere_in_cone (1,:) cell; % ConvexPolyhedron;
  end % End of properties block

  methods
    function this = FlowTransitionGainDigraph(conical_partition, flow_set_cone_ndxs, flow_map_matrix)
      % Create a flow transition gain graph for one or more conical partitions. 
      % 
      % Usage:
      %   
      %   
      %
      % conical_partition_1 = ConicalPartition.fromNumSlices(4);
      % conical_partition_2 = ConicalPartition.fromNumSlices(8);
      % jump_graph = JumpTransitionGainDigraph(conical_partition, jump_specs);
      %   
      %   
      %   
      %   

      arguments(Input, Repeating)
        conical_partition ConicalPartition;
        flow_set_cone_ndxs (1, :) {pwintz.validators.mustBeIndexVector};
        flow_map_matrix (:, :)    {pwintz.validators.mustBeSquare};
      end % End of Input arguments block
      this@TransitionGainDigraph(conical_partition{:});

      % ⋘──────── Validate inputs ────────⋙
      for q = this.mode_indices
        pwintz.assertions.assertAllAreMembers(flow_set_cone_ndxs{q}, conical_partition{q}.cone_indices);
      end

      % ⋘──────── Assign inputs ────────⋙
      this.flow_set_cone_ndxs = flow_set_cone_ndxs;
      this.flow_map_matrix    = flow_map_matrix;

      % ╭────────────────────────────────────────────────╮
      % │             Add edges to the graph             │
      % ╰────────────────────────────────────────────────╯
      for mode_ndx = 1:numel(conical_partition)
        cp = conical_partition{mode_ndx};
        reach_analyzer = FlowReachabilityAnalyzer(cp, flow_set_cone_ndxs{mode_ndx}, flow_map_matrix{q});

        % ⋘──────── Assign Vertex-to-Cone edges ────────⋙
        for vert_ndx = cp.vertex_indices
          for cone_ndx = reach_analyzer.conesDirectlyReachableFromVertex(vert_ndx) 
            [min_gain, max_gain] = reach_analyzer.gainsFromVertexToCone(vert_ndx, cone_ndx);
            this.addEdgeFromVertexToCone(mode_ndx, vert_ndx, mode_ndx, cone_ndx, min_gain, max_gain);
          end
        end

        % ⋘──────── Assign Cone-to-Vertex edges ────────⋙
        for cone_ndx = cp.cone_indices
          for vert_ndx = reach_analyzer.verticesDirectlyReachableFromCone(cone_ndx) 
            [min_gain, max_gain] = reach_analyzer.gainsFromConeToVertex(cone_ndx, vert_ndx);
            this.addEdgeFromConeToVertex(mode_ndx, cone_ndx, mode_ndx, vert_ndx, min_gain, max_gain);
          end
        end

        % ⋘──────── Assign Vertex-to-Vertex edges ────────⋙
        for start_vertex_ndx = cp.ray_indices
          for end_vertex_ndx = reach_analyzer.verticesDirectlyReachableFromVertex(start_vertex_ndx)
            [min_gain, max_gain] = reach_analyzer.gainsFromVertexToVertex(start_vertex_ndx, end_vertex_ndx);
            this.addEdgeFromVertexToVertex(mode_ndx, start_vertex_ndx, mode_ndx, end_vertex_ndx, min_gain, max_gain);

            % Check values.
            last_edge = this.gains_digraph.Edges(end, :);
            last_edge_nodes = this.gains_digraph.Nodes(last_edge.EndNodes, :);
            pwintz.assertions.assertEqual(last_edge_nodes.ModeIndex(1), mode_ndx, leftName="Edge start mode");
            pwintz.assertions.assertEqual(last_edge_nodes.ModeIndex(2), mode_ndx, leftName="Edge end mode");
          end
        end
      end % End of for loop over modes

      % ╭──────────────────────────────────────────╮
      % │             Check properties             │
      % ╰──────────────────────────────────────────╯
      % for n_row = size(this.gains_digraph.Edges, 1)
      %   edge = this.gains_digraph.Edges(n_row, :);
      %   start_node = edge.EndNodes(1);
      %   end_node   = edge.EndNodes(2);
      %   edge_end_ndxs = edge.EndNodes;
      %   nodes = this.gains_digraph.Nodes(edge.EndNodes, :);
      %   this.gains_digraph.Nodes.ModeIndex(edge.EndNodes, :);
      % end
    end % End of constructor

%     % ╭────────────────────────────────────────╮
%     % │             Reachable Sets             │
%     % ╰────────────────────────────────────────╯
%     function [reachable_set_in_cone, reachable_set] = computeReachableSet(this, cone_ndx, P0)
%       % Compute the cone that is reachable from P0 when flowing according to \dot x \in A*C_i, where C_i is given cone .
%       % The "derivative cone" is represented using vertices on the unit sphere, so we scale it up so that we get a set the sufficiently represents "P0 + A*C_i", at least up to a large radius.
%       arguments(Input)
%         this;
%         cone_ndx {pwintz.validators.mustBeIndexScalar};
%         P0 Polytope;
%       end % End of Input arguments block
%       
%       C_i  = this.conical_partition.getCone(cone_ndx);
%       AC_i = this.flow_map_matrix * C_i;
% 
%       pwintz.assertions.assertAll(C_i.containsPoints(P0.vertices), "We have only implemented this function for the case where P0 is a subset of the cone.");
%       
%       % pwintz.strings.format("P0: %D", P0);
% 
%       % ⋘──────── Check for a degenerate case, w ────────⋙
%       % If all of the vertices of P0 are in the boudnary of the cone, and AC_i points "out" of C_i, then the intersection is not full dimensional. Thus, it is hard to numerically find the intersection. 
%       % In this case, the intersection will just be the portion of P0 that is in C_i. 
%       do_any_point_inward = false;
%       for v = this.flow_map_matrix * P0.vertices % For each direction
%         for p = P0.vertices
%           if C_i.atXDoesVPointInward(p, v)
%             msg =pwintz.strings.format("At vertex p = %f, does the direction v = %f points into cone %d? %D", p, v, cone_ndx, C_i.atXDoesVPointInward(p, v));
%             do_any_point_inward = true;
%             break
%           end
%         end
%       end
% 
%       if ~do_any_point_inward
%         reachable_set = P0;
%         return;
%       end
% 
%       reachable_set         = P0 + AC_i;
%       % reachable_set_as_polytope = Polytope.fromConvexHull([P0.vertices + 100*AC_i.rays]);
%       
%       % % # Plot 
%       % figure(1);
%       % clf();
%       % xlim(1.2*[-1, 1]);
%       % ylim(1.2*[-1, 1]);
%       % axis square;
%       % hold on;
%       % AC_i_as_polytope = AC_i.asPolytope();
%       % C_i.plot("FaceColor", "red", "DisplayName", "C_i");
%       % AC_i_as_polytope.plot("FaceColor", "blue", "DisplayName", "A*C_i");
%       % P0.plot("FaceColor", "black", "DisplayName", "P0");
%       % legend();
% 
%       try
%         reachable_set_in_cone = intersection(C_i, reachable_set);
%       catch err
%         % If we can't compute the intersection, that probably means that flows point out of the cone, so the only points in the reachable set are the initial points.
%         reachable_set_in_cone = P0;
%       end
%       
%       % pwintz.strings.format("do_any_point_inward = %d", do_any_point_inward)
%       % if do_any_point_inward
%       %   % reachable_set_in_cone = intersection(reachable_set, C_i);
%       % else
%       %   reachable_set = P0;
%       % end
%     end

    function out = contractEdgesBetweenJumpGraphNodes(this, jump_graph)
      arguments(Input)
        this;
        jump_graph JumpTransitionGainDigraph;
      end % End of Input arguments block
      
    end % End of function

    function transition_graph = contractEdgesToBetweenCones(this, start_cone_ndxs, end_cone_ndxs)
      % Usage
      % `
      % `   new_transition_graph = this.contractEdgesToBetweenCones(start_cone_ndxs_1, end_cone_ndxs_2, start_cone_ndxs_1, end_cone_ndxs_2);
      % `
      arguments(Input)
        this FlowTransitionGainDigraph;
      end % End of Input arguments block
      arguments(Input, Repeating)
        start_cone_ndxs (1, :) double {pwintz.validators.mustBeIndexVector};
        end_cone_ndxs   (1, :) double {pwintz.validators.mustBeIndexVector};
      end % End of Input arguments block
      
      arguments(Output)
        transition_graph (1, 1) TransitionGainDigraph;
      end % End of Output arguments block

      % # Psuedocode:
      % `
      % `    new_graph = new TransitionGainGraph. 
      % `    for each start index and each end index
      % `       search for a path from the start cone to the end cone in the form 
      % `               [cone] -> [boundary] -> ... -> [boundary] -> [cone].
      % `       for each path found, 
      % `          compute the max and min gains
      % `          add edge to new_graph
      % `       end
      % `    end
      % `   

      % Create the underlying transition graph.
      transition_graph = TransitionGainDigraph(this.conical_partitions{:});

      % Check that the right number of repeating arguments are given. 
      pwintz.assertions.assertNumElements(start_cone_ndxs, this.n_modes);
      pwintz.assertions.assertNumElements(end_cone_ndxs,   this.n_modes);

      for mode_ndx = this.mode_indices
        mode_conical_partition = this.conical_partitions{mode_ndx};

        flow_start_cone_ndxs = intersect(start_cone_ndxs{mode_ndx}, this.flow_set_cone_ndxs{mode_ndx});
        flow_end_cone_ndxs   = intersect(  end_cone_ndxs{mode_ndx}, this.flow_set_cone_ndxs{mode_ndx});
        
        % Check that the cone indices are OK.
        pwintz.assertions.assertAllAreMembers(flow_start_cone_ndxs, mode_conical_partition.cone_indices);
        pwintz.assertions.assertAllAreMembers(flow_end_cone_ndxs,   mode_conical_partition.cone_indices);
        
        % Construct the subgraph from that only includes edges between vertices.
        edges_from_verts_to_verts = this.getEdgesFromVerticesToVertices(mode_ndx, mode_ndx);
        vertex_subgraph = digraph(edges_from_verts_to_verts, this.gains_digraph.Nodes);
        
        % edges_not_from_verts_to_verts = [this.getEdgesFromConesToVertices(mode_ndx, mode_ndx); this.getEdgesFromVerticesToCones(mode_ndx, mode_ndx)];
        % faces_subgraph = digraph(edges_not_from_verts_to_verts, this.gains_digraph.Nodes);
        
        % % We will also use edges from vertices to cones...
        % edges_from_verts_to_cones = this.getEdgesFromVerticesToCones();
        % % ...and from cones to verts.
        % edges_from_cones_to_verts = this.getEdgesFromConesToVertices();
        
        % Get all of the vertices that are next to the start and end cones.
        vertices_adjacent_to_start_cones = mode_conical_partition.getVerticesAdjacentToCone(flow_start_cone_ndxs);
        vertices_adjacent_to_end_cones   = mode_conical_partition.getVerticesAdjacentToCone(flow_end_cone_ndxs);
        
        for start_vertex_ndx = vertices_adjacent_to_start_cones
          start_node_ndx = this.nodeIndexOfModeAndVertex(mode_ndx, start_vertex_ndx);
          for end_vertex_ndx = vertices_adjacent_to_end_cones
            end_node_ndx = this.nodeIndexOfModeAndVertex(mode_ndx, end_vertex_ndx);
            
            % fprintf('Finding all paths from vertex %3d to vertex %3d.\n', start_vertex_ndx, end_vertex_ndx);
            % !! Previously we skipped cases where the vertices are equal, but this actually an important case, since we want to include paths like 
            % `   [cone] -> [vertex] -> [cone]
            % !! which only has path lengths of length zero. 
            [paths_nodes, paths_edges] = vertex_subgraph.allpaths(start_node_ndx, end_node_ndx);
            if isempty(paths_nodes)
              % No paths found.
              fprintf('Did not find any paths from vertex %3d (node %d) to vertex %3d (node %d).\n', start_vertex_ndx, start_node_ndx, end_vertex_ndx, end_node_ndx);
              pwintz.assertions.assertNotEqual(start_vertex_ndx, end_vertex_ndx);
              continue
            end
        
            for i_path = pwintz.arrays.cellVectorIndices(paths_nodes)
              path_nodes = paths_nodes{i_path};
              path_edges = paths_edges{i_path};
              
              pwintz.assertions.assertEqual(path_nodes(1),    start_node_ndx, leftName = "first node in path"); %, 'path_nodes=%s must start at start_node_ndx=%d', mat2str(path_nodes), start_node_ndx);
              pwintz.assertions.assertEqual(path_nodes(end),  end_node_ndx, leftName = "last nod in path"); %, 'path_nodes=%s must end at end_node_ndx=%d', mat2str(path_nodes), end_node_ndx);
              % fprintf('Found path for start_cone_ndx=%d, end_cone_ndx=%d, start_vertex_ndx=%d, end_vertex_ndx=%d, i_path=%d\n', start_cone_ndx, end_cone_ndx, start_vertex_ndx, end_vertex_ndx, i_path);
              % fprintf('Found path from vertex %3d (node %d) to vertex %3d (node %d): %s.\n', start_vertex_ndx, start_node_ndx, end_vertex_ndx, end_node_ndx, mat2str(path_nodes));
        
              % Now that we have a path through the vertices, we need to add the start and end edges that connect to the adjacent cones.
              for start_cone_ndx = flow_start_cone_ndxs
                % start_cone_node_ndx = this.nodeIndexOfModeAndCone(mode_ndx, start_cone_ndx);
                start_edge = this.getEdgesFromConeToVertex(mode_ndx, start_cone_ndx, mode_ndx, start_vertex_ndx);
                if isempty(start_edge) 
                  % fprintf('\tNo start_edge from start_cone_ndx=%d (node %d) to start_vertex_ndx=%d. Continuing to next start cone.\n', start_cone_ndx, start_cone_node_ndx, start_vertex_ndx);
                  continue
                else
                  % fprintf('\tFound a start edge from start_cone_ndx=%d (node %d) to start_vertex_ndx=%d. Continuing to next start cone.\n', start_cone_ndx, start_cone_node_ndx, start_vertex_ndx);
                end
                for end_cone_ndx = flow_end_cone_ndxs
                  % end_cone_node_ndx   = this.nodeIndexOfModeAndCone(mode_ndx, end_cone_ndx);
                  end_edge = this.getEdgesFromVertexToCone(mode_ndx, end_vertex_ndx, mode_ndx, end_cone_ndx);
                  if isempty(end_edge) 
                    % fprintf('\t\tNo end_edge from end_vertex_ndx=%d to end_cone_ndx=%d (node %d). Continuing to next end cone.\n', end_vertex_ndx, end_cone_ndx, end_cone_node_ndx);
                    continue
                  end
                  % fprintf('\t\tFound a complete path! edges from end_vertex_ndx=%d to end_cone_ndx=%d (node %d). Continuing to next end cone.\n', end_vertex_ndx, end_cone_ndx, end_cone_node_ndx);
        
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
                  transition_graph.addEdgeFromConeToCone(mode_ndx, start_cone_ndx, mode_ndx, end_cone_ndx, min_gain, max_gain);

                  pwintz.strings.format("Added edge to transition graph (Mode %d, Cone %d) -> (Mode %d, Cone %d) with min_gain = %.2g, max_gain = %.2g.", mode_ndx, start_cone_ndx, mode_ndx, end_cone_ndx, min_gain, max_gain)
                end
              end
            end % End of "path" for loop.
          end % End of for loop
        end % End of for loop
      end % End of "mode_ndx" for loop

    end % End of function



%   end % End of methods block.
% 
%   methods(Access = {?matlab.unittest.TestCase}) % Private block, accessible by unit tests.
  end % End private methods block

end % End of class