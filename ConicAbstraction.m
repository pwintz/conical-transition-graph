classdef ConicAbstraction < handle
  
  %%% TODO: It seems like this class should be merged into ConicalPartition.
  %%% X Option 1: Merge classes, naming the resulting class "ConicaAbstraction".
  %%% * Option 2: Move all of the reachability analysis from ConicalPartion to here. Making ConicalPartion a "data" class.
  %%%  * Details: The ConicalPartition class currently handles constructing two separate partitions: a "derivative" partition and a "state" partition. We could store these in two separate partitions, with the connection being handled by the ConicAbstraction. Choosing this option would allow the ConicalPartition class to be nicely focused on the geometry (which regions are where, how are they connected) without worrying about the dynamics.
  %%% * Option 3: Keep logic for computing local (single-step) flow and jump reachability into ConicalPartition. Then, in ConicalAbstraction, handle the composition of the reachability into graphs for analysis.
  %%%  * Details: This seems like an inferior options. It's distrubtes the handling of the dynamics between multiple classes.
  
  properties(SetAccess = immutable)
    % Define instance constants.
    conical_partition;
    flow_set_region_ndxs;
    jump_set_region_ndxs;
    jump_set_image_region_ndxs;
    flow_graph; % Digraph where each node is a vertex in the conical partition.
    jump_graph; % Digraph where each node is a region in the conical partition.
  end
  
  methods(Static)
    function test(test_function_name)
      test_name = "TestConicAbstraction";
      if nargin() == 1 % If a test function name is given, append it so only that test is run.
        test_name = test_name + "/" + test_function_name;
      end
      results = runtests(test_name);
      fprintf("%d Passed, %d Failed, %d Incomplete.\n", sum([results.Passed]), sum([results.Failed]), sum([results.Incomplete]));
    end % End of function
  end % End static methods block
  
  methods(Static)
    function conic_abstraction = fromAngles(options)
      arguments(Input)
        options.n (1, 1) int32;
        options.flowMapMatrix (2, 2) double; % "A_c" in \dot x = A_c x.
        options.jumpMapMatrix (2, 2) double; % "A_d" in    x^+ = A_d x.
        options.flowSetAngles (1, 2) double;
        options.jumpSetAngles (1, 2) double;
      end % End of Input arguments block
      
      arguments(Output)
        conic_abstraction ConicAbstraction;
      end % End of Output arguments block
      
      A_c = options.flowMapMatrix;
      A_d= options.jumpMapMatrix;
      
      % Create the baseline derivative vertex angles. We'll add more state angles to this list to account for the flow and jump sets.
      derivative_vertex_angles = pwintz.arrays.range("start", 0, "end", 2*pi, "n_values", options.n, "includeEnd", false);
      
      % Compute the image of the jump set under the jump map.
      jump_set_image_vectors = A_d * pwintz.math.angle2UnitVector(options.jumpSetAngles);
      jump_set_image_angles = pwintz.math.atan2(jump_set_image_vectors);
      
      % Create a list of all the aditional state vertex angles we need to include in the ConicalPartition.
      additional_state_vertex_angles = [options.flowSetAngles, options.jumpSetAngles, jump_set_image_angles];
      
      % Create the ConicalPartition.
      conical_partition = ConicalPartition.fromAngles(A_c, derivative_vertex_angles, additional_state_vertex_angles);
      
      
      % Find the ConicalPartition regions that partition the flow set, jump set, and image of the jump set.
      flow_set_region_ndxs       = conical_partition.getRegionsBetweenVerticesFromAngles(options.flowSetAngles(1), options.flowSetAngles(2));
      jump_set_region_ndxs       = conical_partition.getRegionsBetweenVerticesFromAngles(options.jumpSetAngles(1), options.jumpSetAngles(2));
      jump_set_image_region_ndxs = conical_partition.getRegionsBetweenVerticesFromAngles(jump_set_image_angles(1), jump_set_image_angles(2));
      
      % Create the conical partition.
      conic_abstraction = ConicAbstraction(conical_partition, flow_set_region_ndxs, jump_set_region_ndxs, jump_set_image_region_ndxs);
    end % End of function
  end % End static methods block
  
  methods
    
    % ╭─────────────────────────────────────────╮
    % │ ╭─────────────────────────────────────╮ │
    % │ │             Constructor             │ │
    % │ ╰─────────────────────────────────────╯ │
    % ╰─────────────────────────────────────────╯
    function this = ConicAbstraction(conical_partition, flow_set_region_ndxs, jump_set_region_ndxs, jump_set_image_region_ndxs)
      % Function signature should be like: myFunction(<positional arguments>, options)
      arguments
        conical_partition ConicalPartition;
        flow_set_region_ndxs       (1, :) double;
        jump_set_region_ndxs       (1, :) double;
        jump_set_image_region_ndxs (1, :) double;
        % flow_set ConvexPolyhedron; % <- This should be a cone.
      end
      
      this.conical_partition = conical_partition;
      this.flow_set_region_ndxs = flow_set_region_ndxs;
      this.jump_set_region_ndxs = jump_set_region_ndxs;
      this.jump_set_image_region_ndxs = jump_set_image_region_ndxs;
      
      % ╭─────────────────────────────────────────╮
      % │             Construct Graph             │
      % ╰─────────────────────────────────────────╯
      node_table = struct2table(struct(...
        "Name",     arrayfun(@vertexIndex2Name, conical_partition.vertex_indices'),...
        "TeXLabel", arrayfun(@vertexIndex2TexLabel, conical_partition.vertex_indices'),...
        "Index",    conical_partition.vertex_indices',...
        "x",        conical_partition.state_vertices(1, :)',...
        "y",        conical_partition.state_vertices(2, :)'...
        )...
        );
      edge_table = struct2table(struct(...
        "EndNodes",          double.empty(0, 2), ...
        "Weight",            double.empty(0, 1), ...
        "WeightIntervalStr", string.empty(0, 1), ...
        "Length",            double.empty(0, 1) ...
        ) ...
        );
      flow_graph = digraph(edge_table, node_table);
      
      for v0_ndx = conical_partition.vertex_indices
        v0 = conical_partition.getStateSpaceVertex(v0_ndx);
        v0_name = vertexIndex2Name(v0_ndx);
        
        % ⋘────────── Get the neighbors of v0 ──────────⋙
        for v_nb_ndx = conical_partition.getNeighborVertexIndices(v0_ndx)
          
          v_nb = conical_partition.getStateSpaceVertex(v_nb_ndx);
          v_nb_name = vertexIndex2Name(v_nb_ndx);
          % fprintf('Checking if there is an arrow from %s to %s\n', v0_name, v_nb_name);
          conjoining_regions_ndxs = conical_partition.getConjoiningRegionsIndices(v0_ndx, v_nb_ndx);
          assert(numel(conjoining_regions_ndxs) == 1, "We expect exactly one region conjoining a pair of neighboring vertices (in 2D).");
          
          if ~ismember(conjoining_regions_ndxs, this.flow_set_region_ndxs)
            % Region is not in the flow set, so we skip it.
            continue
          end
          
          % ⋘────────── Check if v_nb is reachable from v0 ──────────⋙
          C = conical_partition.getStateCone(conjoining_regions_ndxs);
          AC = conical_partition.getDerivativeCone(conjoining_regions_ndxs);
          % Get the region reachable from v0 when flowing in the conjoing region C under the \dot x = D = AC, which will include any solution to \dot x = Ax in the region.
          R = (v0 + 1e4 * AC) | (1e4 * C);
          
          % ⋘────────── If v_nb is reachable from v0, construct edge──────────⋙
          if ~isempty(R)
            % Get the set of points in ray(v_nb) that is reachable from v0.
            reachable_edge = R.removeVertex(v0).vertices;
            
            % % Check "reachable_edge" is acutally a line segment aligned with ray(v0)
            % TODO: Debug why these assertions fail.
            % pwintz.assertions.assertEqual(size(reachable_edge, 2), 2);
            % pwintz.assertions.assertAlmostEqual(v0 / norm(v0), reachable_edge(:, 1) / norm(reachable_edge(:, 1)));
            % pwintz.assertions.assertAlmostEqual(v0 / norm(v0), reachable_edge(:, 2) / norm(reachable_edge(:, 2)));
            
            w_min = min(vecnorm(reachable_edge));
            w_max = max(vecnorm(reachable_edge));
            assert(w_max >= w_min);
            % interval_weight = [w_min, w_max];
            fprintf('Adding an arrow from %6s to %6s.\tw_min = %s,\tw_max = %s\n', v0_name, v_nb_name, ctg.utils.num2strNear1(w_min), ctg.utils.num2strNear1(w_max));
            
            interval_label = sprintf("[%.2g, %.2g]", w_min, w_max);
            edge_table = table(...
              [v0_name, v_nb_name], ...
              w_max, ...
              interval_label, ...
              norm(v0 - v_nb), ...
              'VariableNames', {'EndNodes', 'Weight', 'WeightIntervalStr', 'Length'});
            flow_graph = flow_graph.addedge(edge_table);
            % graph = graph.addedge(v0_name, v_nb_name, w_max);
          end % End "if ~isempty(R)" block
          % disp("Intersection (R_" + v0v_prev_ndxs + " ∩ v_" + v_prev_ndx + "): ")
          % disp(R_prev.intersectRay(v_prev))
        end % End of for block
        this.flow_graph = flow_graph;
      end % End of for block
      
      region_node_table = struct2table(struct(...
        "Name",     arrayfun(@regionIndex2Name, conical_partition.region_indices'),...
        "TeXLabel", arrayfun(@regionIndex2TexLabel, conical_partition.vertex_indices'),...
        "Index",    conical_partition.vertex_indices',...
        "region",        conical_partition.state_vertices(1, :)',...
        "y",        conical_partition.state_vertices(2, :)'...
        )...
        );
      edge_table = struct2table(struct(...
        "EndNodes",          double.empty(0, 2), ...
        "Weight",            double.empty(0, 1), ...
        "WeightIntervalStr", string.empty(0, 1), ...
        "Length",            double.empty(0, 1) ...
        ) ...
        );
      jump_graph = digraph(edge_table, node_table);
      % for
      %   
      end % End of function
      
      function [v_ndxs, reach_sets_in_cone] = computeDirectlyBackwardReachableVerticesFromSet(this, cone)
        % Find all of the vertices v such that "cone" can be reached from v by flowing through through a single region of the partition.
        % Outputs:
        % * v_ndxs: An array of the vertex indices, arranged horizontally.
        % * reach_sets: An array of the vertex indices, arranged horizontally.
        arguments
          this ConicAbstraction;
          cone ConvexPolyhedron; % <-  We should use a "cone" data type.
        end
        
        % ⋘────────── Find all of the regions that intersect "cone" ──────────⋙
        
      end % End of function
  end % End methods block.
end
  
function v_name = vertexIndex2Name(v_ndx)
  v_name =  "v_" + v_ndx;
end

function label = vertexIndex2TexLabel(v_ndx)
label =  "$v_{" + v_ndx + "}$";
end

function R_name = regionIndex2Name(R_ndx)
R_name =  "R_" + R_ndx;
end

function label = regionIndex2TexLabel(R_ndx)
label =  "$R_{" + R_ndx + "}$";
end