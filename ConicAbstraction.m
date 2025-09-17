classdef ConicAbstraction < handle
  
  properties(SetAccess = immutable)
    % Define instance constants.
    conical_partition; % Partition of state space.
    flow_map_matrix; % "A_c"
    jump_map_matrix; % "A_d"

    flow_set_cone_ndxs         (1, :) int32;
    jump_set_cone_ndxs         (1, :) int32;
    jump_set_image_cone_ndxs   (1, :) int32;
    

    is_cone_in_flow_set        (1, :) logical;
    is_cone_in_jump_set        (1, :) logical;
    is_cone_in_jump_set_image  (1, :) logical;

    n_flow_set_cones           (1, 1);
    n_jump_set_cones           (1, 1);
    n_jump_set_image_cones     (1, 1);

    % can_flow_from_vertex_into_cone is a logical array with dimensions (n_vertices)x(n_cones).
    % The (i,j) entry of can_flow_from_vertex_into_cone is 1 if and only if at vertex i, the flow of \dot x = A_c x does not point out of the cone j, indicating that a flow from v_i can travel into cone_j (or, possibly, travel along the border).
    can_flow_from_vertex_into_cone (:, :); % Logical (n_vertices)x(n_cones)
    can_flow_from_cone_to_vertex   (:, :); % Logical (n_cones)x(n_vertices)

    restricted_reachable_sets_from_unit_sphere_in_cone  (1, :) cell; % Contains ConvexPolyhedrons.             
    reachable_sets_from_unit_sphere_in_cone             (1, :) cell; % Contains ConvexPolyhedrons.  
    directly_reachable_sets_from_vertices               (:, :) cell; % Contains ConvexPolyhedrons.

    flow_transition_graphs (1, :) FlowTransitionGainDigraph;
    jump_transition_graph            TransitionGainDigraph;
    contracted_flow_transition_graph TransitionGainDigraph;
    ctg                              TransitionGainDigraph;

    is_origin_asymptotically_stable logical; % (1, 1) TrueFalseIndederminate;
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
      conical_partition = ConicalPartition.fromUVSphere3D(...
        nLinesOfLongitude=options.nLinesOfLongitude, ...
        nLinesOfLatitude=options.nLinesOfLatitude...
      );
      

      jump_set_center = [-1; 0; 0];
      jump_set_indicator_fnc = @(x) vecnorm(x - jump_set_center, inf) < 0.3;
      flow_set_indicator_fnc = @(x) x(3, :) >= -0.1; % Flow set contains upper halfspace.

      cone_middle_vectors = conical_partition.getConeMiddleVector();

      pwintz.strings.format("ConicAbstraction: Finding cone indices for D in preliminary ConicalPartition.");
      jump_set_cone_ndxs = find(jump_set_indicator_fnc(cone_middle_vectors));

      jump_set_verts_ndxs = conical_partition.getVerticesAdjacentToCone(jump_set_cone_ndxs);
      jump_set_verts = conical_partition.getVertex(jump_set_verts_ndxs);
      jump_set_image_verts =  A_d * jump_set_verts;

      disp("Constructing a new ConicalPartition that includes vertices for G(D).");
      new_vertices = [conical_partition.rays, jump_set_image_verts];
      conical_partition = ConicalPartition(new_vertices);
      cone_middle_vectors = conical_partition.getConeMiddleVector();

      % Find the ConicalPartition cones that partition the flow set, jump set, and image of the jump set.
      pwintz.strings.format("ConicAbstraction: Finding cone indices for flow set C.");
      flow_set_cone_ndxs = find(flow_set_indicator_fnc(cone_middle_vectors));

      pwintz.strings.format("ConicAbstraction: Finding cone indices for D.");
      jump_set_cone_ndxs = find(jump_set_indicator_fnc(cone_middle_vectors));
      % jump_set_middle_vectors = cone_middle_vectors(:, jump_set_cone_ndxs);
      % jump_set_verts_ndxs = conical_partition.getVerticesAdjacentToCone(jump_set_cone_ndxs)

      % jump_set_image_middle_vectors = A_d * jump_set_middle_vectors
      % jump_set_image_cone_ndxs = conical_partition.getConesContainingPoint(cone_middle_vectors(:, jump_set_cone_ndxs))

      is_jump_set_image_cone_ndx = false(1, conical_partition.n_cones);
      for jump_cone_ndx = jump_set_cone_ndxs
        disp("Finding intersections for jump set cone " + jump_cone_ndx);
        cone = conical_partition.getCone(jump_cone_ndx);
        
        [image_cone_ndxs, does_cone_intersect] = conical_partition.getConeIntersections(cone);
        % does_cone_intersect
        is_jump_set_image_cone_ndx(does_cone_intersect) = true;
      end

      % jump_set_verts = conical_partition.getVertex(jump_set_verts_ndxs)
      % jump_set_image_verts =  
    
      % A_d * jump_set_verts

      pwintz.strings.format("ConicAbstraction: Finding cone indices for G(D).");
      jump_set_image_cone_ndxs = find(is_jump_set_image_cone_ndx);

      pwintz.strings.format("ConicAbstraction: Finished finding cone indices for C, D, and G(D).");
      
      % Create the conical abstraction.
      conic_abstraction = ConicAbstraction(conical_partition, flow_set_cone_ndxs, jump_set_cone_ndxs, jump_set_image_cone_ndxs, A_c, A_d);
    end % End of function

    function conic_abstraction = fromAngles2D(options)
      arguments(Input)
        options.maxStateConeAngle      (1, 1) double = 2*pi / 10;
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

      % options.flowSetAngles;
      % options.jumpSetAngles;
      % options_cell = namedargs2cell(options);

      state_cone_angles = generateConeAngles2D(...
        flowSetStartAngle      = options.flowSetAngles(1),  ...
        flowSetEndAngle        = options.flowSetAngles(2),  ...
        jumpSetStartAngle      = options.jumpSetAngles(1),  ...
        jumpSetEndAngle        = options.jumpSetAngles(2),  ...
        maxStateConeAngle      = options.maxStateConeAngle, ...
        maxDerivativeConeAngle = options.maxDerivativeConeAngle);
      
      % Create the ConicalPartition.
      
      pwintz.strings.format("ConicAbstraction: Constructing ConicalPartition from angles.");
      conical_partition = ConicalPartition.fromAngles2D(state_cone_angles);
      
      % Find the ConicalPartition cones that partition the flow set, jump set, and image of the jump set.
      
      jump_set_image_angles = mapAnglesByLinearMap(options.jumpMapMatrix, options.jumpSetAngles);

      pwintz.strings.format("ConicAbstraction: Finding cone indices for C.");
      flow_set_cone_ndxs       = conical_partition.getConesIntersectingArc(options.flowSetAngles(1), options.flowSetAngles(2));
      pwintz.strings.format("ConicAbstraction: Finding cone indices for D.");
      jump_set_cone_ndxs       = conical_partition.getConesIntersectingArc(options.jumpSetAngles(1), options.jumpSetAngles(2));
      pwintz.strings.format("ConicAbstraction: Finding cone indices for G(D).");
      jump_set_image_cone_ndxs = conical_partition.getConesIntersectingArc(jump_set_image_angles(1), jump_set_image_angles(2));
      pwintz.strings.format("ConicAbstraction: Finished finding cone indices for C, D, and G(D).");
      
      % Create the conical abstraction.
      conic_abstraction = ConicAbstraction(conical_partition, flow_set_cone_ndxs, jump_set_cone_ndxs, jump_set_image_cone_ndxs, options.flowMapMatrix, options.jumpMapMatrix);

    end % End of function
  end % End static methods block
  
  
  methods
    % ╭─────────────────────────────────────────╮
    % │ ╭─────────────────────────────────────╮ │
    % │ │             Constructor             │ │
    % │ ╰─────────────────────────────────────╯ │
    % ╰─────────────────────────────────────────╯
    function this = ConicAbstraction(...
        conical_partitions, ...
        flow_specifications, ...
        jump_specifications... 
      )
      arguments
        conical_partitions  (1, :) ConicalPartition {mustBeNonempty};
        flow_specifications (1, :) FlowSpecification {mustBeNonempty};
        % The entry (i,j) specifies a jump from mode i to mode j. (This does not match matrix multiplication, where each each column corresponds to an input component and each row correspondes to an output componenet, but it seems more natural to read jump_specifications(i, j) as the specification of "jumps from i to j" .)
        jump_specifications (:, :) JumpSpecification {pwintz.validators.mustBeSquare};
      end
      %     function this = ConicAbstraction(...
      %         conical_partitions, ...
      %         flow_set_cone_ndxs, ...
      %         jump_set_cone_ndxs, ...
      %         jump_set_image_cone_ndxs, ...
      %         flow_map_matrix, ...
      %         jump_map_matrix  ...
      %       )
      %       arguments
      %         conical_partitions       (1, :) ConicalPartition;
      %         flow_set_cone_ndxs       (1, :) double;
      %         jump_set_cone_ndxs       (1, :) double;
      %         jump_set_image_cone_ndxs (1, :) double;
      % 
      %         % Linear maps Ac and Ad that define \dot x = Ax * x and x^+ = A_d * x.
      %         flow_map_matrix          (:, :) double;
      %         jump_map_matrix          (:, :) double;
      %       end

      n_modes = numel(conical_partitions);
      pwintz.assertions.assertSameSize(conical_partitions, flow_specifications);
      pwintz.assertions.assertSize(jump_specifications, [n_modes, n_modes]);

      % TODO: Check that all of the conical_partitions have the same dimension.
      % dimension = conical_partitions(1).dimension;

      % pwintz.assertions.assertSize(flow_map_matrix, [dimension, dimension]);
      % pwintz.assertions.assertSize(jump_map_matrix, [dimension, dimension]);
      
%       this.conical_partition  = conical_partition;
%       this.flow_set_cone_ndxs = flow_set_cone_ndxs;
%       this.jump_set_cone_ndxs = jump_set_cone_ndxs;
%       this.jump_set_image_cone_ndxs = jump_set_image_cone_ndxs;
% 
%       cone_indices = conical_partition.cone_indices;
%       this.is_cone_in_flow_set       = ismember(cone_indices, flow_set_cone_ndxs);
%       this.is_cone_in_jump_set       = ismember(cone_indices, jump_set_cone_ndxs);
%       this.is_cone_in_jump_set_image = ismember(cone_indices, jump_set_image_cone_ndxs);
% 
%       this.n_flow_set_cones       = numel(this.flow_set_cone_ndxs);
%       this.n_jump_set_cones       = numel(this.jump_set_cone_ndxs);
%       this.n_jump_set_image_cones = numel(this.jump_set_image_cone_ndxs);
% 
%       % ⋘──────── Checks ────────⋙
%       if ~any(this.is_cone_in_jump_set | this.is_cone_in_flow_set)
%         warning("The jump set does not intersect the flow set.");
%       end
%       
%       if ~any(this.is_cone_in_jump_set_image | this.is_cone_in_flow_set)
%         warning("The jump set image does not intersect the flow set.");
%       end
% 
%       % Linear maps Ac and Ad that define \dot x = Ax * x and x^+ = A_d * x.
%       this.flow_map_matrix = flow_map_matrix;
%       this.jump_map_matrix = jump_map_matrix;
% 
% TODO: Restore checks
%       % ⋘──────── Check the indices are all valid ────────⋙
%       % Check that the set of indices for the flow set, jump set, and jump set images are subsets of "conical_partition.cone_indices".
%       pwintz.assertions.assertAllAreMembers(flow_set_cone_ndxs,       cone_indices);
%       pwintz.assertions.assertAllAreMembers(jump_set_cone_ndxs,       cone_indices);
%       pwintz.assertions.assertAllAreMembers(jump_set_image_cone_ndxs, cone_indices);
%       pwintz.assertions.assertUnique(flow_set_cone_ndxs);
%       pwintz.assertions.assertUnique(jump_set_cone_ndxs);
%       pwintz.assertions.assertUnique(jump_set_image_cone_ndxs);

      % ╭────────────────────────────────────────────────────────────────╮
      % │             Find all jump set indices in each mode             │
      % ╰────────────────────────────────────────────────────────────────╯
      pwintz.strings.format("ConicAbstraction: Constructing JumpTransitionGainDigraph.");
      is_jump_set_ndx = cell(1, n_modes);
      for start_node_ndx = 1:n_modes
        start_conical_partition = conical_partition(start_node_ndx);
        is_jump_set_ndx = false(1, start_conical_partition.n_cones);

        for destination_mode_ndx = 1:n_modes
          % The jumps out of mode "node_ndex are given in the "node_ndx" row of jump_specifications.
          jump_spec = jump_specifications(start_node_ndx, destination_mode_ndx);
          pwintz.assertions.assertAllAreMembers(jump_spec.jump_cone_set_ndxs, conical_partitions.cone_indices);
          jump_transition_graph = JumpTransitionGainDigraph(conical_partition, this.jump_set_cone_ndxs, this.jump_map_matrix);

        end

      end


      % ╭──────────────────────────────────────────────╮
      % │             Construct Flow Graphs             │
      % ╰──────────────────────────────────────────────╯
      % * To check that the origin is stable for \dot x = A_c*x, x \in C, we need to check the following:
      %   * For each cycle in flow_graph, the weight is <= 1.
      %   * For each flow set cone,  


      pwintz.strings.format("ConicAbstraction: Constructing FlowTransitionGainDigraphs for each mode.");
      flow_transition_graphs = cell(1, n_modes);
      for mode_ndx = 1:n_modes
        % TODO: Add some sanity checks here (does conical_partition match the flow spec? is dimension OK?)
        flow_spec = flow_specifications(mode_ndx);
        flow_transition_graphs{mode_ndx} = FlowTransitionGainDigraph(conical_partitions(mode_ndx), flow_spec.flow_set_cone_ndxs, flow_spec.flow_map_matrix);

        pwintz.strings.format("ConicAbstraction: Contracting Edges in FlowTransitionGainDigraph.");
        this.contracted_flow_transition_graphs{mode_ndx} = flow_transition_graphs{mode_ndx}.contractEdgesToBetweenCones(this.jump_set_image_cone_ndxs, this.jump_set_cone_ndxs);
      end

      % Assign to property.
      this.flow_transition_graphs = [flow_transition_graphs{:}];


      
      
      % ⋘──────── Construct Jump Graph for Reachability between Cones ────────⋙
      
      % ⋘──────── Construct Conical Transition Graph ────────⋙
      this.ctg = TransitionGainDigraph.union(this.jump_transition_graph, this.contracted_flow_transition_graph);
      pwintz.strings.format("ConicAbstraction: Finding all cycles in CTG.");
      [cycle_min_gains, cycle_max_gains, cycles_nodes, cycles_edges] =  this.ctg.getCycleGains();


      this.is_origin_asymptotically_stable = all(cycle_max_gains < 1) && this.hasStableFlows();

      % for v0_ndx = conical_partition.vertex_indices
      %   v0 = conical_partition.getVertex(v0_ndx);
      %   v0_name = vertexIndex2Name(v0_ndx);
      %   
      %   % ⋘────────── Get the neighbors of v0 ──────────⋙
      %   % !!! Using adjacent nodes works in 2D, but in higher dimensions, to find all of the vertices that are on the boundary of the same cone, we need to use a different approach.
      %   for v_nb_ndx = conical_partition.getVerticesAdjacentToVertex(v0_ndx)
      %     assert(all(v0_ndx ~= v_nb_ndx));
      %     
      %     v_nb      = conical_partition.getVertex(v_nb_ndx);
      %     v_nb_name = vertexIndex2Name(v_nb_ndx);
      %     % fprintf('Checking if there is an arrow from %s to %s\n', v0_name, v_nb_name);
      %     conjoining_cones_ndxs = conical_partition.getStateConeConjoiningTwoVertices(v0_ndx, v_nb_ndx);
      %     if v0_ndx ~= conical_partition.origin_index && v_nb_ndx ~= conical_partition.origin_index
      %       assert(numel(conjoining_cones_ndxs) == 1, "We expect exactly one cone conjoining the pair of vertices v0_ndx = %s, v_nb_ndx = %s (in 2D). Instead there were %d. They are %s.", mat2str(v0_ndx), mat2str(v_nb_ndx), numel(conjoining_cones_ndxs), mat2str(conjoining_cones_ndxs));
      %     end
      %     if ~ismember(conjoining_cones_ndxs, this.flow_set_cone_ndxs)
      %       % Cone is not in the flow set, so we skip it.
      %       continue
      %     end
      %     
      %     % ⋘────────── Check if v_nb is reachable from v0 ──────────⋙
      %     state_cone = this.getStateCone(conjoining_cones_ndxs);
      %     derivative_cone = this.flow_map_matrix * state_cone;
      %     % Get the cone reachable from v0 when flowing in the conjoing cone state_cone under the \dot x = D = derivative_cone, which will include any solution to \dot x = Ax in the cone.
      %     reach_set = (v0 + 1e4 * derivative_cone) | (1e4 * state_cone);
      %     
      %     % ⋘────────── If v_nb is reachable from v0, construct edge──────────⋙
      %     if ~isempty(reach_set)
      %       % Get the set of points in ray(v_nb) that is reachable from v0.
      %       reachable_edge = reach_set.removeVertex(v0).vertices;
      %       
      %       % % Check "reachable_edge" is acutally a line segment aligned with ray(v0)
      %       % TODO: Debug why these assertions fail.
      %       % pwintz.assertions.assertEqual(size(reachable_edge, 2), 2);
      %       % pwintz.assertions.assertAlmostEqual(v0 / norm(v0), reachable_edge(:, 1) / norm(reachable_edge(:, 1)));
      %       % pwintz.assertions.assertAlmostEqual(v0 / norm(v0), reachable_edge(:, 2) / norm(reachable_edge(:, 2)));
      %       
      %       w_min = min(vecnorm(reachable_edge));
      %       w_max = max(vecnorm(reachable_edge));
      %       assert(w_max >= w_min);
      %       % interval_weight = [w_min, w_max];
      %       fprintf('Adding an arrow from %6s to %6s.\tw_min = %s,\tw_max = %s\n', v0_name, v_nb_name, ctg.utils.num2strNear1(w_min), ctg.utils.num2strNear1(w_max));
      %     end % End "if ~isempty(reach_set)" block
      %   end % End of for block
      % end % End of for block
      
      disp("Finished constructing ConicAbstraction");
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
      state_cone = this.conical_partition.getCone(cone_index);
    end % End of function

    function set = getReachableSetInCone(this, cone_index, initial_set)
      
    end % End of function


    % ╭────────────────────────────────────────╮
    % │             Reachable Sets             │
    % ╰────────────────────────────────────────╯
    function [reachable_set_in_cone, reachable_set] = computeReachableSet(this, cone_ndx, P0)
      % Compute the cone that is reachable from P0 when flowing according to \dot x \in A*C_i in a given cone C_i.
      % The "derivative cone" is represented using vertices on the unit sphere, so we scale it up so that we get a set the sufficiently represents "P0 + A*C_i", at least up to a large radius.
      error("Deprecated. See FlowTransitionGainDigraph.computeReachableSet");

      C_i  = this.conical_partition.getCone(cone_ndx);
      AC_i = this.flow_map_matrix * C_i;
      reachable_set         = P0 + AC_i;


      
      reachable_set_in_cone = intersection(reachable_set, C_i);
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
    

    % function reach_set = getReachableSetFromVertex(this, vertex_ndx, min_radius, max_radius, depth)
    %   this.flow_transition_graph.get
    %   reachable_convex_polyhedrons = {};
    %   % if depth = 
    % end
    % 
    % function reach_set = getReachableSetFromCone(this, vertex_ndx, depth)
    %   reachable_convex_polyhedrons = {};
    %   % if depth = 
    % end
    
    % ╭────────────────────────────────────────╮
    % │  ╭──────────────────────────────────╮  │
    % │  │             Plotting             │  │
    % │  ╰──────────────────────────────────╯  │
    % ╰────────────────────────────────────────╯
    function plotVertices(this)
      for i_vertex_ndx = this.conical_partition.ray_indices
        i_vertex = this.conical_partition.getVertex(i_vertex_ndx);
        % Plot the vertex.
        pwintz.plots.plotVector2(i_vertex, plotArgs={'k', "ShowArrowHead", false, "HandleVisibility", "off"});
      end
    end

    function plotCones(this)
      alpha = 0.5;
      % ⋘──────── For each cone ────────⋙
      for cone_ndx = this.conical_partition.cone_indices
          
        % ⋘──────── Plot flow_set ────────⋙
        if ismember(cone_ndx, this.flow_set_cone_ndxs)
          this.conical_partition.plotCone(cone_ndx, "FaceColor", "blue", "FaceAlpha", alpha);
        end

        % ⋘──────── Plot jump_set ────────⋙
        if ismember(cone_ndx, this.jump_set_cone_ndxs)
          this.conical_partition.plotCone(cone_ndx, "FaceColor", "red", "FaceAlpha", alpha);
          % K = this.getStateCone(cone_ndx);
          % G_of_K = this.jump_map_matrix * K;
        end

        % ⋘──────── Plot jump_set_image ────────⋙
        if ismember(cone_ndx, this.jump_set_image_cone_ndxs)
          % disp("Plot jump set cone index #" + cone_ndx);
          this.conical_partition.plotCone(cone_ndx, "FaceColor", [0.9, 0.5, 0], "FaceAlpha", alpha);
        end
      end
    end

    function plotConesReachableFromVertices(this)
      for i_ray_ndx = this.conical_partition.ray_indices
        i_ray = this.conical_partition.getRay(i_ray_ndx);
    
        
        for reachable_cone_ndx = find(this.flow_transition_graph.can_flow_from_ray_into_cone(i_ray_ndx, :))
          % ⋘──────── Plot arrow to middle of reachable cone ────────⋙
          cone_middle_vector = this.conical_partition.getConeMiddleVector(reachable_cone_ndx);
            pwintz.plots.plotVector2(i_ray, cone_middle_vector - i_ray, ...
               plotArgs={"Color", [0, 0.7, 0], "HandleVisibility", "off", "LineWidth", 3, "MaxHeadSize", 3});
        end
      end
    end % End of function

    function plotVerticesReachableFromCones(this)
      for i_cone_ndx = this.conical_partition.cone_indices
        cone_middle_vector = this.conical_partition.getConeMiddleVector(i_cone_ndx);
    
        % ⋘──────── For each eachable vertex, plot an arrow from the cone cone to the vertex ────────⋙
        for reachable_ray_ndx = find(this.flow_transition_graph. can_flow_from_cone_to_ray(i_cone_ndx, :))
          reachable_vertex = this.conical_partition.getVertex(reachable_ray_ndx);
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



% function y = f(options)
%   arguments(Input)
%     options.x
%   end % End of Input arguments block
%    y = g(options);
% end
% 
% function y = g(options)
%   arguments(Input)
%     options.x;
%   end % End of Input arguments block
%   y = options.x^2
% end