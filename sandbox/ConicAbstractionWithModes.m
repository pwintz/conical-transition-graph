classdef ConicAbstractionWithModes < handle
  
  properties(SetAccess = immutable)
    % Define instance constants.

    conic_abstractions;

    jump_maps (:, :) cell;
    jump_set_cones_ndxs         (:, :) cell;
    jump_set_image_cones_ndxs   (:, :) cell;

    n_modes = 2;
%     is_cone_in_flow_set        (1, :) logical;
%     is_cone_in_jump_set        (1, :) logical;
%     is_cone_in_jump_set_image  (1, :) logical;
% 
%     n_flow_set_cones           (1, 1);
%     n_jump_set_cones           (1, 1);
%     n_jump_set_image_cones     (1, 1);

    % flow_graph; % Digraph where each node is a vertex in the conical partition.
    % jump_graph; % Digraph where each node is a cone in the conical partition.

    % can_flow_from_vertex_into_cone is a logical array with dimensions (n_vertices)x(n_cones).
    % The (i,j) entry of can_flow_from_vertex_into_cone is 1 if and only if at vertex i, the flow of \dot x = A_c x does not point out of the cone j, indicating that a flow from v_i can travel into cone_j (or, possibly, travel along the border).
    % can_flow_from_vertex_into_cone (:, :); % Logical (n_vertices)x(n_cones)
    % can_flow_from_cone_to_vertex   (:, :); % Logical (n_cones)x(n_vertices)
    % reachable_cones_from_vertex cell;

%     restricted_reachable_sets_from_unit_sphere_in_cone  (1, :) cell; % Contains ConvexPolyhedrons.             
%     reachable_sets_from_unit_sphere_in_cone             (1, :) cell; % Contains ConvexPolyhedrons.  
%     directly_reachable_sets_from_vertices               (:, :) cell; % Contains ConvexPolyhedrons.
% 
%     flow_transition_graph            FlowTransitionGainDigraph;
%     jump_transition_graph            TransitionGainDigraph;
%     contracted_flow_transition_graph TransitionGainDigraph;
    ctg                              TransitionGainDigraph;

    % is_origin_asymptotically_stable logical; % (1, 1) TrueFalseIndederminate;
  end
  
  methods(Static)
    
    function conic_abstraction = fromUVSphere3D(options)
      arguments(Input)
        options.flowMapMatrix (3, 3) double = rand(3, 3); % "A_c" in \dot x = A_c x.
        options.jumpMapMatrix (3, 3) double = rand(3, 3); % "A_d" in    x^+ = A_d x.
        options.nLinesOfLongitude = 10;
        options.nLinesOfLatitude  = 10;
        options.verbose logical = false;
      end % End of Input arguments block
      
      arguments(Output)
        conic_abstraction ConicAbstraction;
      end % End of Output arguments block
      
      A_c = options.flowMapMatrix;
      A_d = options.jumpMapMatrix;
      ConicAbstraction.checkMatrices(A_c, A_d);

      % Create the ConicalPartition.
      pwintz.strings.format("ConicAbstraction: Constructing preliminary ConicalPartition from UV Sphere (does not include vertices for G(D)).");
      conical_partitions = ConicalPartition.fromUVSphere3D(...
        nLinesOfLongitude=options.nLinesOfLongitude, ...
        nLinesOfLatitude=options.nLinesOfLatitude...
      );
      

      jump_set_center = [-1; 0; 0];
      jump_set_indicator_fnc = @(x) vecnorm(x - jump_set_center, inf) < 0.3;
      flow_set_indicator_fnc = @(x) x(3, :) >= -0.1; % Flow set contains upper halfspace.

      cone_middle_vectors = conical_partitions.getConeMiddleVector();

      pwintz.strings.format("ConicAbstraction: Finding cone indices for D in preliminary ConicalPartition.");
      jump_set_cone_ndxs = find(jump_set_indicator_fnc(cone_middle_vectors));

      jump_set_verts_ndxs = conical_partitions.getVerticesAdjacentToCone(jump_set_cone_ndxs);
      jump_set_verts = conical_partitions.getVertex(jump_set_verts_ndxs);
      jump_set_image_verts =  A_d * jump_set_verts;

      disp("Constructing a new ConicalPartition that includes vertices for G(D).");
      new_vertices = [conical_partitions.rays, jump_set_image_verts];
      conical_partitions = ConicalPartition(new_vertices);
      cone_middle_vectors = conical_partitions.getConeMiddleVector();

      % Find the ConicalPartition cones that partition the flow set, jump set, and image of the jump set.
      pwintz.strings.format("ConicAbstraction: Finding cone indices for flow set C.");
      flow_set_cone_ndxs = find(flow_set_indicator_fnc(cone_middle_vectors));

      pwintz.strings.format("ConicAbstraction: Finding cone indices for D.");
      jump_set_cone_ndxs = find(jump_set_indicator_fnc(cone_middle_vectors))
      % jump_set_middle_vectors = cone_middle_vectors(:, jump_set_cone_ndxs);
      % jump_set_verts_ndxs = conical_partitions.getVerticesAdjacentToCone(jump_set_cone_ndxs)

      % jump_set_image_middle_vectors = A_d * jump_set_middle_vectors
      % jump_set_image_cone_ndxs = conical_partitions.getConesContainingPoint(cone_middle_vectors(:, jump_set_cone_ndxs))

      is_jump_set_image_cone_ndx = false(1, conical_partitions.n_cones);
      for jump_cone_ndx = jump_set_cone_ndxs
        disp("Finding intersections for jump set cone " + jump_cone_ndx);
        cone = conical_partitions.getCone(jump_cone_ndx)
        
        [image_cone_ndxs, does_cone_intersect] = conical_partitions.getConeIntersections(cone)
        does_cone_intersect
        is_jump_set_image_cone_ndx(does_cone_intersect) = true
      end

      % jump_set_verts = conical_partitions.getVertex(jump_set_verts_ndxs)
      % jump_set_image_verts =  
    
      % A_d * jump_set_verts

      pwintz.strings.format("ConicAbstraction: Finding cone indices for G(D).");
      jump_set_image_cone_ndxs = find(is_jump_set_image_cone_ndx);

      pwintz.strings.format("ConicAbstraction: Finished finding cone indices for C, D, and G(D).");
      
      % Create the conical abstraction.
      conic_abstraction = ConicAbstraction(conical_partitions, flow_set_cone_ndxs, jump_set_cone_ndxs, jump_set_image_cone_ndxs, A_c, A_d);
    end % End of function


    function conic_abstraction = fromAngles2D(options)
      arguments(Input)
        options.maxStateConeAngle (1, 1) double = 2*pi / 10;
        options.maxDerivativeConeAngle (1, 1) double = 2*pi / 20;
        options.flowMapMatrix (2, 2) double; % "A_c" in \dot x = A_c x.
        options.jumpMapMatrix (2, 2) double; % "A_d" in    x^+ = A_d x.
        options.flowSetAngles (1, 2) double;
        options.jumpSetAngles (1, 2) double;
        options.verbose logical = false;
      end % End of Input arguments block
      
      arguments(Output)
        conic_abstraction ConicAbstraction;
      end % End of Output arguments block
      
      A_c = options.flowMapMatrix;
      A_d = options.jumpMapMatrix;
      ConicAbstraction.checkMatrices(A_c, A_d);

      ConicAbstraction.checkSetAngles(options.flowSetAngles);
      ConicAbstraction.checkSetAngles(options.jumpSetAngles);

      % Compute the image of the jump set under the jump map.
      jump_set_image_vectors = A_d * pwintz.math.angle2UnitVector(options.jumpSetAngles);
      jump_set_image_angles  = pwintz.math.atan2(jump_set_image_vectors);

      % Sort the two angles in jump_set_image_angles so that the angle from the first to the second is less than pi in the CCW direction. 
      if pwintz.math.angleDiffCCW(jump_set_image_angles, index=1) > pi
        jump_set_image_angles = circshift(jump_set_image_angles, 1);
      end
      assert(pwintz.math.angleDiffCCW(jump_set_image_angles, index=1) <= pi, ...
        'The jump_set_image_angles=%s should be ordered now such that the CCW distance from the first to second entry is < pi.')
      
      % Create the baseline derivative vertex angles. We'll add more state angles to this list to account for the flow and jump sets.
      n_min_derivative_cones = ceil(2*pi / options.maxDerivativeConeAngle);
      derivative_cone_angles = pwintz.arrays.range("start", 0, "end", 2*pi, "n_values", n_min_derivative_cones, "includeEnd", false);
      derivative_cone_vertices = pwintz.math.angle2UnitVector(derivative_cone_angles);
      state_cone_vertices = A_c \ derivative_cone_vertices;
      state_cone_angles = pwintz.math.atan2(state_cone_vertices);

      % Add the angles from the flow set, jump set, and image of the jump.
      state_cone_angles = [state_cone_angles, options.flowSetAngles, options.jumpSetAngles, jump_set_image_angles];
      % Update range to be [0, 2*pi)
      state_cone_angles = mod(state_cone_angles, 2*pi);
      % Sort and remove duplicates.
      state_cone_angles = unique(state_cone_angles, 'sorted');

      % ⋘──────── If any angles are too large, insert needed nodes ────────⋙
      angle_diff = pwintz.math.angleDiffCCW(state_cone_angles);
      extra_angles = [];
      for i = find(angle_diff > options.maxStateConeAngle)
        % We use a crude method for inserting more state cones. We simple take the last angle before the offending interval and step forward at the maximum allowed angle step. 
        extra_angles = [extra_angles, state_cone_angles(i) + (0:options.maxStateConeAngle:angle_diff(i))]; %#ok<AGROW>
      end
      state_cone_angles = unique([state_cone_angles,extra_angles], 'sorted');
      
      % Create the ConicalPartition.
      
      pwintz.strings.format("ConicAbstraction: Constructing ConicalPartition from angles.");
      conical_partitions = ConicalPartition.fromAngles2D(state_cone_angles);
      
      % Find the ConicalPartition cones that partition the flow set, jump set, and image of the jump set.
      
      pwintz.strings.format("ConicAbstraction: Finding cone indices for C.");
      flow_set_cone_ndxs       = conical_partitions.getConesIntersectingArc(options.flowSetAngles(1), options.flowSetAngles(2));
      pwintz.strings.format("ConicAbstraction: Finding cone indices for D.");
      jump_set_cone_ndxs       = conical_partitions.getConesIntersectingArc(options.jumpSetAngles(1), options.jumpSetAngles(2));
      pwintz.strings.format("ConicAbstraction: Finding cone indices for G(D).");
      jump_set_image_cone_ndxs = conical_partitions.getConesIntersectingArc(jump_set_image_angles(1), jump_set_image_angles(2));
      pwintz.strings.format("ConicAbstraction: Finished finding cone indices for C, D, and G(D).");
      
      
      % Create the conical abstraction.
      conic_abstraction = ConicAbstraction(conical_partitions, flow_set_cone_ndxs, jump_set_cone_ndxs, jump_set_image_cone_ndxs, A_c, A_d);

    end % End of function
  end % End static methods block
  
  
  methods
    % ╭─────────────────────────────────────────╮
    % │ ╭─────────────────────────────────────╮ │
    % │ │             Constructor             │ │
    % │ ╰─────────────────────────────────────╯ │
    % ╰─────────────────────────────────────────╯
    function this = ConicAbstractionWithModes(conical_partitions, flow_set_cone_ndxs, jump_set_cone_ndxs, jump_set_image_cone_ndxs, flow_map_matrix, jump_map_matrix)
      % Function signature should be like: myFunction(<positional arguments>, options)
      arguments
        conical_partitions ConicalPartition;
        flow_set_cone_ndxs       (1, :) double;
        jump_set_cone_ndxs       (1, :) double;
        jump_set_image_cone_ndxs (1, :) double;

        % Linear maps Ac and Ad that define \dot x = Ax * x and x^+ = A_d * x.
        flow_map_matrix          (:, :) double;
        jump_map_matrix          (:, :) double;
      end

      dimension = conical_partitions.dimension;

      pwintz.assertions.assertSize(flow_map_matrix, [dimension, dimension]);
      pwintz.assertions.assertSize(jump_map_matrix, [dimension, dimension]);
      
      this.conical_partitions  = conical_partitions;
      this.flow_set_cone_ndxs = flow_set_cone_ndxs;
      this.jump_set_cone_ndxs = jump_set_cone_ndxs;
      this.jump_set_image_cone_ndxs = jump_set_image_cone_ndxs;

      cone_indices = conical_partitions.cone_indices;
      this.is_cone_in_flow_set       = ismember(cone_indices, flow_set_cone_ndxs);
      this.is_cone_in_jump_set       = ismember(cone_indices, jump_set_cone_ndxs);
      this.is_cone_in_jump_set_image = ismember(cone_indices, jump_set_image_cone_ndxs);

      this.n_flow_set_cones       = numel(this.flow_set_cone_ndxs);
      this.n_jump_set_cones       = numel(this.jump_set_cone_ndxs);
      this.n_jump_set_image_cones = numel(this.jump_set_image_cone_ndxs);

      % ⋘──────── Checks ────────⋙
      if ~any(this.is_cone_in_jump_set | this.is_cone_in_flow_set)
        warning("The jump set does not intersect the flow set.");
      end
      
      if ~any(this.is_cone_in_jump_set_image | this.is_cone_in_flow_set)
        warning("The jump set image does not intersect the flow set.");
      end

      % Linear maps Ac and Ad that define \dot x = Ax * x and x^+ = A_d * x.
      this.flow_map_matrix = flow_map_matrix;
      this.jump_map_matrix = jump_map_matrix;

      % ⋘──────── Check the indices are all valid ────────⋙
      % Check that the set of indices for the flow set, jump set, and jump set images are subsets of "conical_partitions.cone_indices".
      pwintz.assertions.assertAllAreMembers(flow_set_cone_ndxs,       cone_indices);
      pwintz.assertions.assertAllAreMembers(jump_set_cone_ndxs,       cone_indices);
      pwintz.assertions.assertAllAreMembers(jump_set_image_cone_ndxs, cone_indices);
      pwintz.assertions.assertUnique(flow_set_cone_ndxs);
      pwintz.assertions.assertUnique(jump_set_cone_ndxs);
      pwintz.assertions.assertUnique(jump_set_image_cone_ndxs);

      % ╭──────────────────────────────────────────────╮
      % │             Construct Flow Graph             │
      % ╰──────────────────────────────────────────────╯
      % * To check that the origin is stable for \dot x = A_c*x, x \in C, we need to check the following:
      %   * For each cycle in flow_graph, the weight is <= 1.
      %   * For each flow set cone,  
      pwintz.strings.format("ConicAbstraction: Constructing FlowTransitionGainDigraph.");
      this.flow_transition_graph = FlowTransitionGainDigraph(conical_partitions, flow_set_cone_ndxs, flow_map_matrix);

      pwintz.strings.format("ConicAbstraction: Contracting Edges in FlowTransitionGainDigraph.");
      this.contracted_flow_transition_graph = this.flow_transition_graph.contractEdgesToBetweenCones(this.jump_set_image_cone_ndxs, this.jump_set_cone_ndxs);
      
      % ⋘──────── Construct Jump Graph for Reachability between Cones ────────⋙
      pwintz.strings.format("ConicAbstraction: Constructing JumpTransitionGainDigraph.");
      
      this.jump_transition_graph = JumpTransitionGainDigraph(conical_partitions, this.jump_set_cone_ndxs, this.jump_set_image_cone_ndxs, this.jump_map_matrix);
      
      % ⋘──────── Construct Conical Transition Graph ────────⋙
      this.ctg = TransitionGainDigraph.union(this.jump_transition_graph, this.contracted_flow_transition_graph);
      pwintz.strings.format("ConicAbstraction: Finding all cycles in CTG.");
      [cycle_min_gains, cycle_max_gains, cycles_nodes, cycles_edges] =  this.ctg.getCycleGains();


      this.is_origin_asymptotically_stable = all(cycle_max_gains < 1) && this.hasStableFlows();


      disp("Returning early before the end of ConicAbstraction constructor.");
      return
      for v0_ndx = conical_partitions.vertex_indices
        v0 = conical_partitions.getVertex(v0_ndx);
        v0_name = vertexIndex2Name(v0_ndx);
        
        % ⋘────────── Get the neighbors of v0 ──────────⋙
        % !!! Using adjacent nodes works in 2D, but in higher dimensions, to find all of the vertices that are on the boundary of the same cone, we need to use a different approach.
        for v_nb_ndx = conical_partitions.getVerticesAdjacentToVertex(v0_ndx)
          assert(all(v0_ndx ~= v_nb_ndx));
          
          v_nb      = conical_partitions.getVertex(v_nb_ndx);
          v_nb_name = vertexIndex2Name(v_nb_ndx);
          % fprintf('Checking if there is an arrow from %s to %s\n', v0_name, v_nb_name);
          conjoining_cones_ndxs = conical_partitions.getStateConeConjoiningTwoVertices(v0_ndx, v_nb_ndx);
          if v0_ndx ~= conical_partitions.origin_index && v_nb_ndx ~= conical_partitions.origin_index
            assert(numel(conjoining_cones_ndxs) == 1, "We expect exactly one cone conjoining the pair of vertices v0_ndx = %s, v_nb_ndx = %s (in 2D). Instead there were %d. They are %s.", mat2str(v0_ndx), mat2str(v_nb_ndx), numel(conjoining_cones_ndxs), mat2str(conjoining_cones_ndxs));
          end
          if ~ismember(conjoining_cones_ndxs, this.flow_set_cone_ndxs)
            % Cone is not in the flow set, so we skip it.
            continue
          end
          
          % ⋘────────── Check if v_nb is reachable from v0 ──────────⋙
          state_cone = this.getStateCone(conjoining_cones_ndxs);
          derivative_cone = this.flow_map_matrix * state_cone;
          % Get the cone reachable from v0 when flowing in the conjoing cone state_cone under the \dot x = D = derivative_cone, which will include any solution to \dot x = Ax in the cone.
          reach_set = (v0 + 1e4 * derivative_cone) | (1e4 * state_cone);
          
          % ⋘────────── If v_nb is reachable from v0, construct edge──────────⋙
          if ~isempty(reach_set)
            % Get the set of points in ray(v_nb) that is reachable from v0.
            reachable_edge = reach_set.removeVertex(v0).vertices;
            
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
            
            % interval_label = sprintf("[%.2g, %.2g]", w_min, w_max);
            % edge_table = table(...
            %   [v0_name, v_nb_name], ...
            %   w_max, ...
            %   interval_label, ...
            %   norm(v0 - v_nb), ...
            %   'VariableNames', {'EndNodes', 'Weight', 'WeightIntervalStr', 'Length'});
            % flow_graph = flow_graph.addedge(edge_table);
            flow_transition_graph.addEdgeFromVertexToVertex(v0_ndx, v_nb_ndx, w_min, w_max);

            % graph = graph.addedge(v0_name, v_nb_name, w_max);
          end % End "if ~isempty(reach_set)" block
          % disp("Intersection (R_" + v0v_prev_ndxs + " ∩ v_" + v_prev_ndx + "): ")
          % disp(R_prev.intersectRay(v_prev))
        end % End of for block
        % this.flow_graph = flow_graph;
      end % End of for block
      
      % node_table = struct2table(struct(...
      %   "Name",     arrayfun(@coneIndex2Name, conical_partitions.cone_indices'),...
      %   "TeXLabel", arrayfun(@coneIndex2TexLabel, conical_partitions.vertex_indices'),...
      %   "Index",    conical_partitions.vertex_indices',...
      %   "x",        conical_partitions.state_vertices(1, :)',...
      %   "y",        conical_partitions.state_vertices(2, :)'...
      % ));
      % edge_table = struct2table(struct(...
      %   "EndNodes",          double.empty(0, 2), ...
      %   "Weight",            double.empty(0, 1), ...
      %   "WeightIntervalStr", string.empty(0, 1), ...
      %   "Length",            double.empty(0, 1) ...
      % ));
      % jump_graph = digraph(edge_table, node_table);
      % for
      %

      pwintz.strings.format("Finished constructing %s", this)
    end % End of function
    
    function [v_ndxs, reach_sets_in_cone] = computeDirectlyBackwardReachableVerticesFromSet(this, cone)
      % Find all of the vertices v such that "cone" can be reached from v by flowing through through a single cone of the partition.
      % Outputs:
      % * v_ndxs: An array of the vertex indices, arranged horizontally.
      % * reach_sets: An array of the vertex indices, arranged horizontally.
      arguments
        this ConicAbstraction;
        cone ConvexPolyhedron; % <-  We should use a "cone" data type.
      end
      
      % ⋘────────── Find all of the cones that intersect "cone" ──────────⋙
      
    end % End of function

    function state_cone = getStateCone(this, cone_index)
      state_cone = this.conical_partitions.getCone(cone_index);
    end % End of function

    function set = getReachableSetInCone(this, cone_index, initial_set)
      
    end % End of function


    % ╭────────────────────────────────────────╮
    % │             Reachable Sets             │
    % ╰────────────────────────────────────────╯
    function [reachable_set_in_cone, reachable_set] = computeReachableSet(this, cone_ndx, P0)
      % Compute the cone that is reachable from P0 when flowing according to \dot x \in A*C_i, where C_i is given cone .
      % The "derivative cone" is represented using vertices on the unit sphere, so we scale it up so that we get a set the sufficiently represents "P0 + A*C_i", at least up to a large radius.
      C_i  = this.conical_partitions.getCone(cone_ndx);
      AC_i = this.flow_map_matrix * C_i;
      reachable_set         = P0 + 1e5*AC_i;
      reachable_set_in_cone = intersection(reachable_set, 1e5*C_i);
    end

    % ╭────────────────────────────────────────╮
    % │  ╭──────────────────────────────────╮  │
    % │  │             Analysis             │  │
    % │  ╰──────────────────────────────────╯  │
    % ╰────────────────────────────────────────╯
    function has_stable_flows = hasStableFlows(this)
      % arguments(Output)
      %   has_stable_flows TrueFalseIndederminate;
      % end % End of Input arguments block
      arguments(Output)
        has_stable_flows logical; % True if stable, false if indeterminate.
      end % End of Output arguments block
      for flow_set_ndx = this.flow_set_cone_ndxs
        R = this.flow_transition_graph.restricted_reachable_sets_from_unit_sphere_in_cone{flow_set_ndx};
        assert(~isempty(R), "All of the cones in the flow set should have a nonempty reachable set.");
        if max(vecnorm(R.vertices)) > 1e3
          has_stable_flows = false; % TrueFalseIndederminate.indeterminate;
          return
        end
      end
      [~, cycles_edges] = this.flow_transition_graph.getVertexCycles();
      assert(isvector(cycles_edges), "Expected cycles_edges to be 1-dimensional. Instead its size was %s.", mat2str(size(cycles_edges)));
      for i_cycle = 1:numel(cycles_edges)
        cycle_edge_ndxs = cycles_edges{i_cycle};
        edges = this.flow_transition_graph.getEdgeRow(cycle_edge_ndxs);
        cycle_max_gain = prod(edges.MaxGain);
        if cycle_max_gain > 1
          has_stable_flows = false; % TrueFalseIndederminate.indeterminate;
          return
        end
      end
      has_stable_flows = true; % TrueFalseIndederminate.true;
    end % End of function

    % ╭─────────────────────────────────────────────────────╮
    % │  ╭───────────────────────────────────────────────╮  │
    % │  │             Reachability Analysis             │  │
    % │  ╰───────────────────────────────────────────────╯  │
    % ╰─────────────────────────────────────────────────────╯
    

    function reach_set = getReachableSetFromVertex(this, vertex_ndx, min_radius, max_radius, depth)
      this.flow_transition_graph.get
      reachable_convex_polyhedrons = {};
      % if depth = 
    end
    
    function reach_set = getReachableSetFromCone(this, vertex_ndx, depth)
      reachable_convex_polyhedrons = {};
      % if depth = 
    end
    
    % ╭────────────────────────────────────────╮
    % │  ╭──────────────────────────────────╮  │
    % │  │             Plotting             │  │
    % │  ╰──────────────────────────────────╯  │
    % ╰────────────────────────────────────────╯
    function plotVertices(this)
      for i_vertex_ndx = this.conical_partitions.ray_indices
        i_vertex = this.conical_partitions.getVertex(i_vertex_ndx);
        % Plot the vertex.
        pwintz.plots.plotVector2(i_vertex, plotArgs={'k', "ShowArrowHead", false, "HandleVisibility", "off"});
      end
    end

    function plotCones(this)
      alpha = 0.5;
      % ⋘──────── For each cone ────────⋙
      for cone_ndx = this.conical_partitions.cone_indices
          
        % ⋘──────── Plot flow_set ────────⋙
        if ismember(cone_ndx, this.flow_set_cone_ndxs)
          this.conical_partitions.plotCone(cone_ndx, "FaceColor", "blue", "FaceAlpha", alpha);
        end

        % ⋘──────── Plot jump_set ────────⋙
        if ismember(cone_ndx, this.jump_set_cone_ndxs)
          this.conical_partitions.plotCone(cone_ndx, "FaceColor", "red", "FaceAlpha", alpha);
          % K = this.getStateCone(cone_ndx);
          % G_of_K = this.jump_map_matrix * K;
        end

        % ⋘──────── Plot jump_set_image ────────⋙
        if ismember(cone_ndx, this.jump_set_image_cone_ndxs)
          disp("Plot jump set cone index" + cone_ndx);
          this.conical_partitions.plotCone(cone_ndx, "FaceColor", [0.9, 0.5, 0], "FaceAlpha", alpha);
        end
      end
    end

    function plotConesReachableFromVertices(this)
      for i_ray_ndx = this.conical_partitions.ray_indices
        i_ray = this.conical_partitions.getRay(i_ray_ndx);
    
        
        for reachable_cone_ndx = find(this.flow_transition_graph.can_flow_from_ray_into_cone(i_ray_ndx, :))
          % ⋘──────── Plot arrow to middle of reachable cone ────────⋙
          cone_middle_vector = this.conical_partitions.getConeMiddleVector(reachable_cone_ndx);
            pwintz.plots.plotVector2(i_ray, cone_middle_vector - i_ray, ...
               plotArgs={"Color", [0, 0.7, 0], "HandleVisibility", "off", "LineWidth", 3, "MaxHeadSize", 3});
        end
      end
    end % End of function

    function plotVerticesReachableFromCones(this)
      for i_cone_ndx = this.conical_partitions.cone_indices
        cone_middle_vector = this.conical_partitions.getConeMiddleVector(i_cone_ndx);
    
        % ⋘──────── For each eachable vertex, plot an arrow from the cone cone to the vertex ────────⋙
        for reachable_vertex_ndx = find(this.can_flow_from_cone_to_vertex(i_cone_ndx, :))
          reachable_vertex = this.conical_partitions.getVertex(reachable_vertex_ndx);
          pwintz.plots.plotVector2(cone_middle_vector, reachable_vertex - cone_middle_vector, ....
              plotArgs={"Color", [0, 0.2, 0.6], "HandleVisibility", "off", "LineWidth", 3, "MaxHeadSize", 3});
        end
      end
    end % End of function

    function plot_handle= plotFlowGraph(this)
      % p = 
      this.flow_transition_graph.plot();
%       graph = this.flow_graph;
%       p = plot(graph, ...
%         "XData", graph.Nodes.x, ...
%         "YData", graph.Nodes.y, ...
%         "NodeColor", "black", ...
%         "EdgeAlpha", 1, ...
%         "MarkerSize", 8, ...
%         "LineWidth", 10, ...
%         'NodeLabel', graph.Nodes.TeXLabel, ...
%         ...'EdgeLabel', graph.Edges.WeightIntervalStr, ...
%         "NodeFontSize", 14, ...
%         ... "EdgeFontSize", 10 * graph.Edges.Length / max(graph.Edges.Length), ...
%         "Interpreter", "latex"...
%       );
% 
%       geq_1_color = [1 0 0];
%       lt_1_color  = [0 0 1];
%       edge_color  = (graph.Edges.Weight >= 1) * geq_1_color + (graph.Edges.Weight < 1) * lt_1_color;
% 
%       % For the plotted weights, we want to illustrate that the importance aspect of an edges weight is its distance above or below 1. We use abs(log(weight)) to normalize. This way, weights equal to 1/2 and 2 are plotted with the same width.
%       weight = abs(log(graph.Edges.Weight));
%       weight = 6 * weight / max(weight) + 1;
%       p.LineWidth = weight;
%       p.ArrowSize = 8*sqrt(p.LineWidth);
%       p.EdgeColor = edge_color;

      if nargout() == 1
        plot_handle = p;
      end
    end % End of function

  end % End methods block.
  
  % ╭──────────────────────────────────────────────╮
  % │  ╭────────────────────────────────────────╮  │
  % │  │             Value Checking             │  │
  % │  ╰────────────────────────────────────────╯  │
  % ╰──────────────────────────────────────────────╯
  methods(Static, Access = private)
    function checkMatrices(A_c, A_d)
      % Check flow map matrix
      assert(isnumeric(A_c));
      assert(isreal(A_c));
      assert(ismatrix(A_c));
      assert(size(A_c,1) == size(A_c,2), "The matrix A_c must be square");

      % Check jump map matrix
      assert(isnumeric(A_d));
      assert(isreal(A_d));
      assert(ismatrix(A_d));
      assert(size(A_d,1) == size(A_d,2), "The matrix A_d must be square");

      % For the flow map matrix, check that it is invertible
      assert(cond(A_c) < 1e9, "The matrix A_c must be well-conditioned but had a condition number of %8.2g", cond(A_c))
      
      % Check the dimension is OK. 
      assert(ismember(size(A_c, 1), [2, 3]), "Only 2D and 3Dimplemented, currently. A_c=%s", mat2str(A_c));
    end % End of function

    function checkSetAngles(angles)
      input_name = inputname(1);
      pwintz.assertions.assertSize(angles, [1, 2]);
      assert(all(angles >= 0), "Angles used to define the set boundaries of ""%s"" must be nonnegative", input_name);
      assert(all(angles < 2*pi), "Angles used to define the set boundaries of ""%s"" must be less than 2*pi", input_name);
      assert(angles(2) >= angles(1), "The end angle used to define ""%s"" must be larger than the start angle", input_name);
    end % End of function
  end % End static methods block
end

function v_name = vertexIndex2Name(v_ndx)
  v_name = sprintf("v_%d", v_ndx);
end

function label = vertexIndex2TexLabel(v_ndx)
  label = sprintf("$v_{%d}$", v_ndx);
  % label =  "$v_{" + v_ndx + "}$";
end

function R_name = coneIndex2Name(R_ndx)
  R_name =  "R_" + R_ndx;
end

function label = coneIndex2TexLabel(R_ndx)
  label =  "$R_{" + R_ndx + "}$";
end