classdef ConicAbstraction < handle
  
  properties(SetAccess = immutable)
    % Define instance constants.
    conical_partition
    vertex_indices
    graph
  end
  properties
    % Define instance variables.
    
  end
  methods
    
    % Constructor
    function this = ConicAbstraction(conical_partition)
      % Function signature should be like: myFunction(<positional arguments>, options)
      arguments
        conical_partition ConicalPartition
      end
      this.conical_partition = conical_partition;
      this.vertex_indices = 1:conical_partition.n_state_vertices;

      % ╭─────────────────────────────────────────╮
      % │             Construct Graph             │
      % ╰─────────────────────────────────────────╯
      n_nodes = conical_partition.n_state_vertices;
      adjacency_mat = zeros(n_nodes); % Start with no edges.
      node_table = table(...
        arrayfun(@vertexIndex2Name, this.vertex_indices'),...
        arrayfun(@vertexIndex2TexLabel, this.vertex_indices'),...
        conical_partition.state_vertices(1, :)',...
        conical_partition.state_vertices(2, :)',...
        'VariableNames', {'Name', 'TeXLabel', 'x', 'y'}...
      );
      edge_table = table(...
                double.empty(0, 2), ...
                double.empty(0, 1), ...
                string.empty(0, 1), ... 
                double.empty(0, 1), ...
                'VariableNames', {'EndNodes', 'Weight', 'WeightIntervalStr', 'Length'})
      graph = digraph(edge_table, node_table);
      graph.Nodes
      
      for v0_ndx = this.vertex_indices
        % ⋘────────── Get the nerighbors of a vertex ──────────⋙
        v0 = conical_partition.getStateSpaceVertex(v0_ndx);
        v0_name = vertexIndex2Name(v0_ndx);
        % v0_nbd_ndxs = conical_partition.getNeighborVertexIndices(v0_ndx);

        for v_nb_ndx = conical_partition.getNeighborVertexIndices(v0_ndx)

          v_nb = conical_partition.getStateSpaceVertex(v_nb_ndx);
          v_nb_name = vertexIndex2Name(v_nb_ndx);
          fprintf('Checking if there is an arrow from %s to %s\n', v0_name, v_nb_name);
          conjoining_regions_ndxs = conical_partition.getConjoiningRegionsIndices(v0_ndx, v_nb_ndx);
          assert(numel(conjoining_regions_ndxs) == 1, "We expect exactly one region conjoining a pair of neighboring vertices (in 2D).")

          % ⋘────────── Check if v_nb is reachable from v0 ──────────⋙
          D = conical_partition.getDerivativeCone(conjoining_regions_ndxs);
          C = conical_partition.getStateCone(conjoining_regions_ndxs);
          % Get the region reachable from v0 when flowing in the conjoing region C under the \dot x = D = AC, which will include any solution to \dot x = Ax in the region.
          R = (v0 + 1e4 * D) | (1e4 * C);
          if ~isempty(R)
            % Get the set of points in ray(v_nb) that is reachable from v0.
            reachable_edge = R.removeVertex(v0).vertices
            vecnorm(reachable_edge)
            w_min = min(vecnorm(reachable_edge));
            w_max = max(vecnorm(reachable_edge));
            assert(w_max >= w_min)
            % new_edge = table(... 
            %   []
            %   [5 9; 3 6], {'on' 'off'}', [10 20]', ...
            %           'VariableNames',{'EndNodes','Power','Weight'});
            interval_weight = [w_min, w_max];
            fprintf('Adding an arrow from %s to %s. w_min = %.2g, w_max = %.2g\n', v0_name, v_nb_name, w_min, w_max);

            interval_label = sprintf("[%.2g, %.2g]", w_min, w_max);
            edge_table = table(...
                [v0_name, v_nb_name], ...
                w_max, ...
                interval_label, ... 
                norm(v0 - v_nb), ...
                'VariableNames', {'EndNodes', 'Weight', 'WeightIntervalStr', 'Length'});
            graph = graph.addedge(edge_table);
            % graph = graph.addedge(v0_name, v_nb_name, w_max);
          end
          % disp("Intersection (R_" + v0v_prev_ndxs + " ∩ v_" + v_prev_ndx + "): ")
          % disp(R_prev.intersectRay(v_prev))
        end
        this.graph = graph;
      end

    end
  end
end

function v_name = vertexIndex2Name(v_ndx)
 v_name =  "v_" + v_ndx;
end

function label = vertexIndex2TexLabel(v_ndx)
 label =  "$v_{" + v_ndx + "}$";
end