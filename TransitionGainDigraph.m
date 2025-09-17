classdef TransitionGainDigraph < handle
  % A class for storing the transitions between cones and vertices in a conical partition. Each transition has associated gains.
  % ` runtests TestTransitionGainDigraph

  properties(SetAccess=private, GetAccess = protected)
    % Define private variables.
    % !! The gains_digraph is mutable because we add edges to it. 
    % !! The nodes of the graph are fixed, however.
    gains_digraph (1,1) digraph; 
  end % End of private properties block

  properties(SetAccess = immutable, GetAccess = protected)
    % Define private variables.
    conical_partitions (1, :) cell; % cell of ConicalPartition;
    n_modes      (1, 1);
    mode_indices (1, :);
  end % End of private properties block

  properties(SetAccess = immutable, GetAccess = private) % Define private variables.
    % Array that contains the graph node indices. 
    % Entry (i, j) contains the index for mode i and vertices j.
    node_ndxs_of_mode_and_vertices (:,1) cell; 
    % Array that contains the graph node indices. 
    % Entry (i, j) contains the index for mode i and cones j.
    node_ndxs_of_mode_and_cones    (:,1) cell; 

    node_ndxs_of_vertices (:,1) {pwintz.validators.mustBeIndexVector}; 
    node_ndxs_of_cones    (:,1) {pwintz.validators.mustBeIndexVector}; 
  end % End of private properties block

%   methods(Static)
%     function gain_graph = empty(array_size)
%       if nargin() ~= 0
%         error("TransitionGainDigraph:NotImplemented","Creating an empty array of TransitionGainDigraphs is not implemented");
%       end
% 
%       
%     end % End of function
%   end % End static methods block

  methods(Static)
    function union_tgd = union(left, right)
      pwintz.assertions.assertEqual(left.gains_digraph.Nodes, right.gains_digraph.Nodes);
      
      % conical_partitions_cell = pwintz.arrays.mat2cellColumn(left.conical_partitions)
      union_tgd = TransitionGainDigraph(left.conical_partitions{:});
      union_edges = [
        left.gains_digraph.Edges;
        right.gains_digraph.Edges;
      ];
      % union_tgd.gains_digraph.Edges = union_edges;
      union_tgd.gains_digraph = union_tgd.gains_digraph.addedge(union_edges);
    end % End of function
  end % End static methods block

  methods
    
    % Constructor
    function this = TransitionGainDigraph(conical_partitions)
      arguments(Input, Repeating)
        conical_partitions (1, 1) ConicalPartition;
      end % End of Input arguments block
      
      n_modes = numel(conical_partitions);

      % n_cones = conical_partitions.n_cones
      % sum(conical_partitions.n_cones) + sum(conical_partitions.n_vertices)

      node_names                     = cell(n_modes, 1);
      mode_ndxs                      = cell(n_modes, 1);
      vertex_ndxs                    = cell(n_modes, 1);
      cone_ndxs                      = cell(n_modes, 1);
      vertex_tex_labels              = cell(n_modes, 1);
      cone_tex_labels                = cell(n_modes, 1);
      objects                        = cell(n_modes, 1);
      tex_label                      = cell(n_modes, 1);
      plot_position                  = cell(n_modes, 1);

      node_ndxs_of_mode_and_vertices = cell(n_modes, 1);
      node_ndxs_of_mode_and_cones    = cell(n_modes, 1);

      last_node_ndx_of_prior_modes = 0;
      for q_mode = 1:n_modes
        conical_partition = conical_partitions{q_mode};
        mode_str = sprintf("Mode %2d, ", q_mode);
        % vert_str = sprintf("Vertex %2d", conical_partition.vertex_indices');
        % cone_str = sprintf("Mode %2d", q_mode);
        node_names{q_mode} = [mode_str + "Vertex " + conical_partition.vertex_indices'; mode_str + "Cone " + conical_partition.cone_indices'];

        mode_ndxs{q_mode}   = repmat(q_mode, size(node_names{q_mode}));
        vertex_ndxs{q_mode} = [conical_partition.vertex_indices'; 0 * conical_partition.cone_indices'];
        cone_ndxs{q_mode}   = [0*conical_partition.vertex_indices'; conical_partition.cone_indices'];
  
        node_ndxs_of_mode_and_vertices{q_mode} = last_node_ndx_of_prior_modes + find(vertex_ndxs{q_mode});
        node_ndxs_of_mode_and_cones{q_mode}    = last_node_ndx_of_prior_modes + find(cone_ndxs{q_mode});

        % Increase the value of "last_node_ndx_of_prior_modes" by the number of nodes in the current mode, which is the number of cones plus the number of vertices.
        last_node_ndx_of_prior_modes = last_node_ndx_of_prior_modes + conical_partition.n_cones + conical_partition.n_vertices;
  
        vertex_tex_labels{q_mode} = vertexIndicesToTexLabels(q_mode, conical_partition.vertex_indices);
        cone_tex_labels{q_mode}   = coneIndicesToTexLabels(q_mode, conical_partition.cone_indices);
  
        objects{q_mode} = [num2cell(conical_partition.vertices, 1), num2cell(conical_partition.cones)]';
        tex_label{q_mode} = [vertex_tex_labels{q_mode}; cone_tex_labels{q_mode}];

        % ⋘──────── Find the plot positions for each cone and vertices ────────⋙
        switch q_mode
          case 1
            mode_plot_center = [-1.2; 0];
          case 2
            mode_plot_center = [+1.2; 0];
          otherwise
            error("TransitionGainDigraph:UNEXPECTED_CASE", "Unexpected case: q_mode=%s.", q_mode);
        end
        cone_centers = arrayfun(@(cone) mode_plot_center + centroid(cone.rays')', conical_partition.cones, "UniformOutput", false);
        vert_cell_array = num2cell(mode_plot_center + conical_partition.vertices, 1);
        
        plot_position{q_mode} = transpose([vert_cell_array, cone_centers]);
      end

      Names        = vertcat(node_names{:});
      ModeIndex    = vertcat(mode_ndxs{:});
      VertexIndex  = vertcat(vertex_ndxs{:});
      ConeIndex    = vertcat(cone_ndxs{:});
      Objects      = vertcat(objects{:});
      TeXLabel     = vertcat(tex_label{:});
      PlotPosition = vertcat(plot_position{:});
      PlotPosition = cell2mat(PlotPosition')'; % Make PlotPosition into a tall array.

      this.node_ndxs_of_mode_and_vertices = node_ndxs_of_mode_and_vertices;
      this.node_ndxs_of_mode_and_cones    = node_ndxs_of_mode_and_cones;

      node_table = pwintz.tables.makeTable(...
        "Names", Names, ...
        "Objects", Objects, ...
        "PlotPosition", PlotPosition, ...
        "ModeIndex", ModeIndex, ...
        "VertexIndex", VertexIndex, ...
        "ConeIndex", ConeIndex, ...
        "TeXLabel", TeXLabel...
      );

      edge_table = struct2table(struct(...
        "EndNodes", double.empty(0, 2), ...
        "Label", string.empty(0, 1), ...
        "MaxGain",  double.empty(0, 1), ...
        "MinGain",  double.empty(0, 1)  ...
      ));
      this.gains_digraph = digraph(edge_table, node_table);
      this.conical_partitions = conical_partitions;

      this.n_modes = n_modes;
      this.mode_indices = 1:n_modes;

      % ╭──────────────────────────────────────────╮
      % │             Check Properties             │
      % ╰──────────────────────────────────────────╯
      
      % ⋘──────── Check all of the node indices, that they aren't cone indices ────────⋙
      for node_ndx = vertcat(node_ndxs_of_mode_and_vertices{:})'
        node_row = this.getNodeRow(node_ndx);
        pwintz.assertions.assertEqual(node_row.ConeIndex, 0, leftName=pwintz.strings.format("ConeIndex of node %d", node_ndx));
      end
      
      
      % ⋘──────── Check all of the node indices, that they aren't cone indices ────────⋙
      for node_ndx = vertcat(node_ndxs_of_mode_and_cones{:})'
        node_row = this.getNodeRow(node_ndx);
        pwintz.assertions.assertEqual(node_row.VertexIndex, 0, leftName=pwintz.strings.format("VertexIndex of node %d", node_ndx));
      end
      % pwintz.strings.format("node_ndxs_of_mode_and_vertices: %s", vertcat(node_ndxs_of_mode_and_vertices{:})')
      % pwintz.strings.format("   node_ndxs_of_mode_and_cones: %s", vertcat(node_ndxs_of_mode_and_cones{:})')
      % this.gains_digraph.Nodes


      return

%       % ╭────────────────────────────────────────────────╮
%       % │             Construct gain_digraph             │
%       % ╰────────────────────────────────────────────────╯
%       if exist("gain_digraph", "var")
%         this.gains_digraph = gains_digraph;
%         error("Do we use this conditional anywhere?");
%       else
% 
%         % node_table = table(...
%         %   Names,          ...
%         %   Objects,        ...
%         %   PlotPosition,   ...
%         %   VertexIndex,    ...
%         %   ConeIndex,      ...
%         %   TeXLabel        ...
%         % )
%         % node_table
% 
%         edge_table = struct2table(struct(...
%           "EndNodes", double.empty(0, 2), ...
%           ..."Weight",   double.empty(0, 1), ...
%           "MaxGain",  double.empty(0, 1), ...
%           "MinGain",  double.empty(0, 1)  ...
%         ));
%         this.gains_digraph = digraph(edge_table, node_table);
%       end % End of if block.
%       this.conical_partition = conical_partition;
    end % End of constructor.
  end % End methods block

  methods
    % ╭──────────────────────────────────────────────────╮
    % │  ╭────────────────────────────────────────────╮  │
    % │  │             Query Graph Properties         │  │
    % │  ╰────────────────────────────────────────────╯  │
    % ╰──────────────────────────────────────────────────╯
    % ╭─────────────────────────────────────╮
    % │             Graph Sizes             │
    % ╰─────────────────────────────────────╯
    function n = numEdges(this)
      n = size(this.gains_digraph.Edges, 1);
    end % End of function
    function n = numNodes(this)
      n = size(this.gains_digraph.Nodes, 1);
    end % End of function
    function n = numVertexNodes(this, mode_ndx)
      n = numel(this.node_ndxs_of_mode_and_vertices{mode_ndx});
      % n = numel(this.gains_digraph.Edges);
    end % End of function
    function n = numConeNodes(this, mode_ndx)
      n = numel(this.node_ndxs_of_mode_and_cones{mode_ndx});
    end % End of function

    % ╭────────────────────────────────────────────────╮
    % │             Get Items By Reference             │
    % ╰────────────────────────────────────────────────╯
    function vertex_row = getVertexRow(this, mode_ndx, vertex_ndx)
      arguments(Input)
        this;
        mode_ndx;
        vertex_ndx; 
      end % End of Input arguments block
      
      node_ndx = this.nodeIndexOfModeAndVertex(mode_ndx, vertex_ndx);
      vertex_row = this.gains_digraph.Nodes(node_ndx, :);
    end % End of function

    function cone_row = getConeRow(this, mode_ndx, cone_ndx)
      arguments(Input)
        this;
        mode_ndx;
        cone_ndx; 
      end % End of Input arguments block
      
      node_ndx = this.nodeIndexOfModeAndCone(mode_ndx, cone_ndx);
      cone_row = this.gains_digraph.Nodes(node_ndx, :);
    end % End of function

    function [edge_row, start_node, end_node] = getEdgeRow(this, edge_ndx)
      arguments(Input)
        this;
        edge_ndx; 
      end % End of Input arguments block
      
      edge_row = this.gains_digraph.Edges(edge_ndx, :);

      start_node = this.gains_digraph.Nodes(edge_row.EndNodes(1), :);
      end_node   = this.gains_digraph.Nodes(edge_row.EndNodes(2), :);
    end % End of function

    function edges = getEdges(this)
      arguments(Input)
        this;
      end % End of Input arguments block
      
      edges = this.gains_digraph.Edges;
    end % End of function

    function nodes = getNodes(this)
      arguments(Input)
        this;
      end % End of Input arguments block
      
      nodes = this.gains_digraph.Nodes;
    end % End of function


    % ╭──────────────────────────────────────────╮
    % │             Has Edge Methods             │
    % ╰──────────────────────────────────────────╯
    function has_edge = hasEdgeFromVertexToCone(this, start_mode_ndx, start_vertex_ndx, end_mode_ndx, end_cone_ndx)
      arguments(Input)
        this;
        start_mode_ndx;
        start_vertex_ndx; 
        end_mode_ndx;
        end_cone_ndx;
      end % End of Input arguments block
      
      start_node_ndx = this.nodeIndexOfModeAndVertex(start_mode_ndx, start_vertex_ndx);
      end_node_ndx   = this.nodeIndexOfModeAndCone(    end_mode_ndx,     end_cone_ndx);
      has_edge = this.hasEdgeBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function

    function has_edge = hasEdgeFromVertexToVertex(this, start_mode_ndx, start_vertex_ndx, end_mode_ndx, end_vertex_ndx)
      arguments(Input)
        this;
        start_mode_ndx;
        start_vertex_ndx; 
        end_mode_ndx;
        end_vertex_ndx;
      end % End of Input arguments block
      
      start_node_ndx = this.nodeIndexOfModeAndVertex(start_mode_ndx, start_vertex_ndx);
      end_node_ndx = this.nodeIndexOfModeAndVertex(end_mode_ndx, end_vertex_ndx);
      has_edge = this.hasEdgeBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function

    function has_edge = hasEdgeFromConeToCone(this, start_mode_ndx, start_cone_ndx, end_mode_ndx, end_cone_ndx)
      arguments(Input)
        this;
        start_mode_ndx;
        start_cone_ndx; 
        end_mode_ndx;
        end_cone_ndx;
      end % End of Input arguments block
      
      start_node_ndx = this.nodeIndexOfModeAndCone(start_mode_ndx, start_cone_ndx);
      end_node_ndx = this.nodeIndexOfModeAndCone(end_mode_ndx, end_cone_ndx);
      has_edge = this.hasEdgeBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function

    function has_edge = hasEdgeFromConeToVertex(this, start_mode_ndx, start_cone_ndx, end_mode_ndx, end_vertex_ndx)
      arguments(Input)
        this;
        start_mode_ndx;
        start_cone_ndx; 
        end_mode_ndx  ;
        end_vertex_ndx;
      end % End of Input arguments block
      
      start_node_ndx = this.nodeIndexOfModeAndCone(start_mode_ndx, start_cone_ndx);
      end_node_ndx = this.nodeIndexOfModeAndVertex(end_mode_ndx, end_vertex_ndx);
      has_edge = this.hasEdgeBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function

    % ╭─────────────────────────────────────────╮
    % │             Get Graph Edges             │
    % ╰─────────────────────────────────────────╯
    function edge_rows = getEdgesFromVertexToCone(this, start_mode_ndx, start_vertex_ndxs, end_mode_ndx, end_cone_ndxs)
      % Get all of the edges that start at given indices in a (single) given mode and end at given indices in a (single) given mode.
      
      % assert(all(start_vertex_ndxs <= size(this.node_ndxs_of_mode_and_vertices, 1)), "start_vertex_ndxs=%s must be less than %d", mat2str(start_vertex_ndxs), size(this.node_ndxs_of_mode_and_vertices, 1));
      % assert(all(end_cone_ndxs <= size(this.node_ndxs_of_mode_and_cones, 1)), "end_cone_ndxs=%s must be less than %d", mat2str(end_cone_ndxs), size(this.node_ndxs_of_mode_and_cones, 1));
      start_node_ndx = this.nodeIndexOfModeAndVertex(start_mode_ndx, start_vertex_ndxs);
      end_node_ndx   = this.nodeIndexOfModeAndCone(end_mode_ndx, end_cone_ndxs);
      edge_rows      = this.getEdgesBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function
    
    function edge_rows = getEdgesFromConeToVertex(this, start_mode_ndx, start_cone_ndx, end_mode_ndx, end_vertex_ndx)
      start_node_ndx = this.nodeIndexOfModeAndCone(start_mode_ndx, start_cone_ndx);
      end_node_ndx   = this.nodeIndexOfModeAndVertex(end_mode_ndx, end_vertex_ndx);
      edge_rows      = this.getEdgesBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function
    
    function edge_rows = getEdgesFromConeToCone(this, start_mode_ndx, start_cone_ndx, end_mode_ndx, end_cone_ndx)    
      start_node_ndx = this.nodeIndexOfModeAndCone(start_mode_ndx, start_cone_ndx);
      end_node_ndx   = this.nodeIndexOfModeAndCone(end_mode_ndx, end_cone_ndx);  
      edge_rows      = this.getEdgesBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function
    
    function edge_rows = getEdgesFromVertexToVertex(this, start_mode_ndx, start_vertex_ndx, end_mode_ndx, end_vertex_ndx)
      start_node_ndx = this.nodeIndexOfModeAndVertex(start_mode_ndx, start_vertex_ndx);
      end_node_ndx   = this.nodeIndexOfModeAndVertex(end_mode_ndx, end_vertex_ndx);
      edge_rows      = this.getEdgesBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function

    % ╭───────────────────────────────────────╮
    % │         All Edges, by Type            │
    % ╰───────────────────────────────────────╯
    function edge_rows = getEdgesFromVerticesToVertices(this, start_mode_ndx, end_mode_ndx)

      if this.n_modes == 1
        start_mode_ndx = 1;
        end_mode_ndx = 1;
      end

      all_edges = this.gains_digraph.Edges;
      edge_start_node_ndxs = all_edges.EndNodes(:, 1);
      edge_end_node_ndxs   = all_edges.EndNodes(:, 2);
      
      does_edge_start_at_vertex = ismember(edge_start_node_ndxs, this.node_ndxs_of_mode_and_vertices{start_mode_ndx});
      does_edge_end_at_vertex   = ismember(edge_end_node_ndxs,   this.node_ndxs_of_mode_and_vertices{end_mode_ndx});

      edge_rows = all_edges(does_edge_start_at_vertex & does_edge_end_at_vertex, :);
      
    end % End of function

    function edge_rows = getEdgesFromVerticesToCones(this, start_mode_ndx, end_mode_ndx)
      if this.n_modes == 1
        start_mode_ndx = 1;
        end_mode_ndx = 1;
      end
      all_edges = this.gains_digraph.Edges;
      edge_start_node_ndxs = all_edges.EndNodes(:, 1);
      edge_end_node_ndxs   = all_edges.EndNodes(:, 2);
      
      does_edge_start_at_vertex = ismember(edge_start_node_ndxs, this.node_ndxs_of_mode_and_vertices{start_mode_ndx});
      does_edge_end_at_cone   = ismember(edge_end_node_ndxs,   this.node_ndxs_of_mode_and_cones{end_mode_ndx});
      
      edge_rows = all_edges(does_edge_start_at_vertex & does_edge_end_at_cone, :);
    end % End of function
    
    function edge_rows = getEdgesFromConesToVertices(this, start_mode_ndx, end_mode_ndx)
      if this.n_modes == 1
        start_mode_ndx = 1;
        end_mode_ndx = 1;
      end
      all_edges = this.gains_digraph.Edges;
      edge_start_node_ndxs = all_edges.EndNodes(:, 1);
      edge_end_node_ndxs   = all_edges.EndNodes(:, 2);
      
      does_edge_start_at_cone = ismember(edge_start_node_ndxs, this.node_ndxs_of_mode_and_cones{start_mode_ndx});
      does_edge_end_at_vertex   = ismember(edge_end_node_ndxs, this.node_ndxs_of_mode_and_vertices{end_mode_ndx});
      
      edge_rows = all_edges(does_edge_start_at_cone & does_edge_end_at_vertex, :);
    end % End of function

    function edge_rows = getEdgesFromConesToCones(this, start_mode_ndx, end_mode_ndx)
      if this.n_modes == 1
        start_mode_ndx = 1;
        end_mode_ndx = 1;
      end
      all_edges = this.gains_digraph.Edges;
      edge_start_node_ndxs = all_edges.EndNodes(:, 1);
      edge_end_node_ndxs   = all_edges.EndNodes(:, 2);
      
      does_edge_start_at_cone = ismember(edge_start_node_ndxs, this.node_ndxs_of_mode_and_cones{start_mode_ndx});
      does_edge_end_at_cone   = ismember(edge_end_node_ndxs,   this.node_ndxs_of_mode_and_cones{end_mode_ndx});
      
      edge_rows = all_edges(does_edge_start_at_cone & does_edge_end_at_cone, :);
      assert(~isempty(edge_rows));
    end % End of function

    % ╭────────────────────────────────────────╮
    % │  ╭──────────────────────────────────╮  │
    % │  │             Analysis             │  │
    % │  ╰──────────────────────────────────╯  │
    % ╰────────────────────────────────────────╯
    % ╭──────────────────────────────────────╮
    % │             Reachability             │
    % ╰──────────────────────────────────────╯
    function vertex_reach_table = getReachableSetFromVertex(this, vertex_ndx, depth, cumulativeMinGain, cumulativeMaxGain, vertex_reach_table, search_depth)
      assert(ismember(nargin(), [2, 7]), "Expected 2 or 7 arguments. Instead had %d", nargin());
      if nargin() == 2
        % if ~exist("vertex_reach_table", "var")
        depth = 0;
        cumulativeMinGain = 1.0;
        cumulativeMaxGain = 1.0;
        vertex_reach_table = struct2table(struct(...
          "vertex", vertex_ndx, ...
          "depth", depth, ...
          "cumulativeMinGain", cumulativeMinGain,  ...
          "cumulativeMaxGain", cumulativeMaxGain  ...
        ));
        search_depth = 20; % How many times to recurse??
      else 
        assert(nargin() == 7, "Expected 1 input or 7. Instead had %d.", nargin());
        if depth == search_depth
          return
        end 
      end % End of if block.
      depth = depth + 1;
      successors = this.gains_digraph.successors(vertex_ndx)';
      if isempty(successors)
        return
      end
      for next_vertex_ndx = this.gains_digraph.successors(vertex_ndx)'
        this.getEdgesFromVertexToVertex(vertex_ndx, next_vertex_ndx);
        
        vertex_reach_table = this.getReachableSetFromVertex(next_vertex_ndx, depth, cumulativeMinGain, cumulativeMaxGain, vertex_reach_table, search_depth);
      end
    end

    % ╭──────────────────────────────────────────╮
    % │             Paths and Cycles             │
    % ╰──────────────────────────────────────────╯
    function [paths_node_ndxs, paths_edge_ndxs, paths_edge_rows] = getPathsBetweenVertices(this, start_vertex_ndx, end_vertex_ndx)
      start_node_ndx = this.nodeIndexOfModeAndVertex(start_vertex_ndx);
      end_node_ndx   = this.nodeIndexOfModeAndVertex(end_vertex_ndx);
      [paths_node_ndxs, paths_edge_ndxs, paths_edge_rows] = this.getPathsBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function

    function [paths_node_ndxs, paths_edge_ndxs, paths_edge_rows] = getPathsBetweenCones(this, start_cone_ndx, end_cone_ndx)
      start_node_ndx = this.node_ndxs_of_mode_and_cones(start_cone_ndx);
      end_node_ndx   = this.node_ndxs_of_mode_and_cones(end_cone_ndx);
      [paths_node_ndxs, paths_edge_ndxs, paths_edge_rows] = this.getPathsBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function

    function [cycles_nodes, cycles_edges] = getCycles(this)
      [cycles_nodes, cycles_edges] = this.gains_digraph.allcycles()
    end % End of function

    function [cycles_nodes, cycles_edges] = getVertexCycles(this)
      vertex_graph = subgraph(this.gains_digraph, this.getAllVertexNodes())  ;
      [cycles_nodes, cycles_edges] = vertex_graph.allcycles();
    end % End of function

    function [cycle_min_gains, cycle_max_gains, cycles_nodes, cycles_edges] = getCycleGains(this)
      [cycles_nodes, cycles_edges] = this.getCycles();
      cycle_min_gains = nan(size(cycles_edges));
      cycle_max_gains = nan(size(cycles_edges));

      assert(isvector(cycles_edges), "Expected cycles_edges to be 1-dimensional. Instead its size was %s.", mat2str(size(cycles_edges)));
      for i_cycle = 1:numel(cycles_edges)
        cycle_edges = cycles_edges{i_cycle};
        
        cycle_path_min_gains = this.gains_digraph.Edges.MinGain(cycle_edges);
        cycle_path_max_gains = this.gains_digraph.Edges.MaxGain(cycle_edges);
        cycle_min_gains(i_cycle) = prod(cycle_path_min_gains);
        cycle_max_gains(i_cycle) = prod(cycle_path_max_gains);
      end
    end % End of function

    function [min_gain, max_gain] = getPathGains(this, path_edge_ndxs)
      arguments(Input)
        this TransitionGainDigraph;
        path_edge_ndxs (:, 1) double {mustBeInteger,mustBePositive};
      end % End of Input arguments block
      
      arguments(Output)
        min_gain (1, 1) double;
        max_gain (1, 1) double;
      end % End of Output arguments block
      
      edge_rows = this.getEdgeRow(path_edge_ndxs);
      min_gains = edge_rows.MinGain;
      max_gains = edge_rows.MaxGain;
      min_gain = prod(min_gains);
      max_gain = prod(max_gains);
      % this.gains_digraph.Edges(path_edge_ndxs)
    end % End of function

    % ╭────────────────────────────────────────╮
    % │  ╭──────────────────────────────────╮  │
    % │  │             Plotting             │  │
    % │  ╰──────────────────────────────────╯  │
    % ╰────────────────────────────────────────╯
    function plot(this, varargin)
      graph = this.gains_digraph;
      positions = graph.Nodes.PlotPosition;
      pwintz.assertions.assertNumColumns(positions, 2);
      blank_labels = repmat("", size(graph.Nodes.TeXLabel));
      p = plot(graph, ...
        "XData", positions(:, 1), ...
        "YData", positions(:, 2), ...
        "NodeColor", "black", ...
        "EdgeAlpha", 0.5, ...
        "Marker", "none", ...
        ..."MarkerSize", 8, ...
        "LineWidth", 1, ...
        ...'NodeLabel', graph.Nodes.Names, ...
        ...'NodeLabel', graph.Nodes.TeXLabel, ...
        'NodeLabel', blank_labels, ...
        'NodeLabelMode', 'manual', ...
        ...'EdgeLabel', graph.Edges.WeightIntervalStr, ...
        "NodeFontSize", 12, ...
        "EdgeColor", "black", ...
        ... "EdgeFontSize", 10 * graph.Edges.Length / max(graph.Edges.Length), ...
        "Interpreter", "latex", ...
        varargin{:}...
      );
      
      % !! Using the TeX labels does not work for unknown reasons. The result is no label at all, so for now we will just use the default "node number" labels.
      ...p.NodeLabel = graph.Nodes.TeXLabel;

      if this.numEdges() > 0
        % ⋘──────── Set color of edges based on whether values are positive or negative ────────⋙
        % geq_1_color = [1 0 0];
        % lt_1_color  = [0 0 1];
        % edge_color  = (graph.Edges.MaxGain >= 1) * geq_1_color + (graph.Edges.MaxGain < 1) * lt_1_color;
        % p.EdgeColor = edge_color;
        
        % ⋘──────── Set the width of the lines based on the weights ────────⋙
        % For the plotted weights, we want to illustrate that the importance aspect of an edges weight is its distance above or below 1. We use abs(log(weight)) to normalize. This way, weights equal to 1/2 and 2 are plotted with the same width.
        weight = abs(log(graph.Edges.MaxGain));
        weight = 6 * (weight + 1) / (max(weight) + 1) + 1;
        % p.LineWidth = weight;

        % ⋘──────── Set arrow size ────────⋙
        p.ArrowSize = 10*sqrt(p.LineWidth);
      end

    end % End of function


    % ╭─────────────────────────────────────────╮
    % │  ╭───────────────────────────────────╮  │
    % │  │             Subgraphs             │  │
    % │  ╰───────────────────────────────────╯  │
    % ╰─────────────────────────────────────────╯
    
%     function vertex_subgraph = getGrapWithOnlyEdgesBetweenVertices(this)
%       arguments(Output)
%         vertex_subgraph (1, 1) digraph;
%       end % End of Output arguments block
%       this.gains_digraph.Edges()
% 
%       ismember(this.gains_digraph.Edges.EndNodes, this.node_ndxs_of_mode_and_cones)
% 
%       this.gains_digraph.rmedge()
%       vertex_subgraph = this.gains_digraph.subgraph(this.node_ndxs_of_mode_and_vertices);
%     end % End of function

%     function vertex_subgraph = getVertexSubgraph(this)
%       % Generate a subgraph that has the same nodes but only has edges between vertex nodes.
%       arguments(Output)
%         vertex_subgraph (1, 1) digraph;
%       end % End of Output arguments block
%       
%       all_nodes     = this.gains_digraph.Nodes;
%       vertex_edges  = this.getEdgesFromVerticesToVertices();
%       vertex_subgraph = digraph(vertex_edges, all_nodes);
%       % vertex_subgraph = this.gains_digraph.subgraph(this.node_ndxs_of_mode_and_vertices);
%     end % End of function
% 
%     function vertex_subgraph = getConeSubgraph(this)
%       arguments(Output)
%         vertex_subgraph (1, 1) digraph;
%       end % End of Output arguments block
%       
%       vertex_subgraph = this.gains_digraph.subgraph(this.node_ndxs_of_mode_and_cones);
%     end % End of function

    % ╭──────────────────────────────────────────────────────────╮
    % │  ╭────────────────────────────────────────────────────╮  │
    % │  │             Overload Builtin Functions             │  │
    % │  ╰────────────────────────────────────────────────────╯  │
    % ╰──────────────────────────────────────────────────────────╯
    function string_rep = char(this)
      % This method defines the string representation that is 
      % inserted when using "%s" in a string format, such as 
      % fprintf("%s", polyhedron) or sprintf("%s", polyhedron); 
      % It does not change the "disp" output. 
      arguments(Output)
        string_rep char; % Ensure output is cast to char, even if you create a string.
      end
      % string_rep = sprintf("%s with %d nodes (%d vertex nodes and %d cone nodes) and %d edges.", class(this), this.numNodes(), this.numVertexNodes(), this.numConeNodes(), this.numEdges());
      string_rep = sprintf("%s with %d modes, %d nodes (?? vertex nodes and ?? cone nodes) and %d edges.", class(this), this.n_modes, this.numNodes(), this.numEdges());
    end

    function disp(this)
      fprintf('%s\n', this);
      pwintz.strings.format("Nodes (%z)", this.gains_digraph.Nodes);
      disp(this.gains_digraph.Nodes);
      pwintz.strings.format("Edges (%z)", this.gains_digraph.Edges);
      disp(this.gains_digraph.Edges);
    end % End of function

    % ╭──────────────────────────────────────────────────────╮
    % │  ╭────────────────────────────────────────────────╮  │
    % │  │             Print Info About Graph             │  │
    % │  ╰────────────────────────────────────────────────╯  │
    % ╰──────────────────────────────────────────────────────╯
    function printInfo(this)
      disp(this);
      disp(this.gains_digraph);
    end % End of function
  end % End of methods block


  methods(Access= {?TransitionGainDigraph, ?matlab.unittest.TestCase}) % Define protected methods.
    % ╭─────────────────────────────────────────╮
    % │  ╭───────────────────────────────────╮  │
    % │  │             Add Edges             │  │
    % │  ╰───────────────────────────────────╯  │
    % ╰─────────────────────────────────────────╯
    function this = addEdgeFromVertexToCone(this, start_mode_ndx, start_vertex_ndx, end_mode_ndx, end_cone_ndx, min_gain, max_gain)
      arguments(Input)
        this;
        start_mode_ndx   (1, 1) {pwintz.validators.mustBeIndexScalar};
        start_vertex_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        end_mode_ndx     (1, 1) {pwintz.validators.mustBeIndexScalar};
        end_cone_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        min_gain (1, 1) double {mustBeNonnegative};
        max_gain (1, 1) double {mustBeNonnegative};
      end % End of Input arguments block
      
      start_node_ndx = this.nodeIndexOfModeAndVertex(start_mode_ndx, start_vertex_ndx);
      end_node_ndx   = this.nodeIndexOfModeAndCone(end_mode_ndx, end_cone_ndx);
      this.addEdgeBetweenNodes(start_node_ndx, end_node_ndx, min_gain, max_gain);
    end % End of function

  
    function this = addEdgeFromConeToVertex(this, start_mode_ndx, start_cone_ndx, end_mode_ndx, end_vertex_ndx, min_gain, max_gain)
      arguments(Input)
        this;
        start_mode_ndx   (1, 1) {pwintz.validators.mustBeIndexScalar};
        start_cone_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        end_mode_ndx     (1, 1) {pwintz.validators.mustBeIndexScalar};
        end_vertex_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        min_gain (1, 1) double {mustBeNonnegative};
        max_gain (1, 1) double {mustBeNonnegative};
      end % End of Input arguments block
      
      start_node_ndx = this.nodeIndexOfModeAndCone(start_mode_ndx, start_cone_ndx);
      end_node_ndx   = this.nodeIndexOfModeAndVertex(end_mode_ndx, end_vertex_ndx);
      this.addEdgeBetweenNodes(start_node_ndx, end_node_ndx, min_gain, max_gain);
    end % End of function
  
    function this = addEdgeFromConeToCone(this, start_mode_ndx, start_cone_ndx, end_mode_ndx, end_cone_ndx, min_gain, max_gain)    
      arguments(Input)
        this;
        start_mode_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        start_cone_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        end_mode_ndx   (1, 1) {pwintz.validators.mustBeIndexScalar};
        end_cone_ndx   (1, 1) {pwintz.validators.mustBeIndexScalar};
        min_gain (1, 1) double {mustBeNonnegative};
        max_gain (1, 1) double {mustBeNonnegative};
      end % End of Input arguments block
      start_node_ndx = this.nodeIndexOfModeAndCone(start_mode_ndx, start_cone_ndx);
      end_node_ndx   = this.nodeIndexOfModeAndCone(end_mode_ndx,   end_cone_ndx);  
      this.addEdgeBetweenNodes(start_node_ndx, end_node_ndx, min_gain, max_gain);
    end % End of function
  
    function this = addEdgeFromVertexToVertex(this, start_mode_ndx, start_vertex_ndx, end_mode_ndx, end_vertex_ndx, min_gain, max_gain)
      arguments(Input)
        this;
        start_mode_ndx   (1, 1) {pwintz.validators.mustBeIndexScalar};
        start_vertex_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        end_mode_ndx     (1, 1) {pwintz.validators.mustBeIndexScalar};
        end_vertex_ndx   (1, 1) {pwintz.validators.mustBeIndexScalar};
        min_gain (1, 1) double {mustBeNonnegative};
        max_gain (1, 1) double {mustBeNonnegative};
      end % End of Input arguments block
      start_node_ndx = this.nodeIndexOfModeAndVertex(start_mode_ndx, start_vertex_ndx);
      end_node_ndx   = this.nodeIndexOfModeAndVertex(end_mode_ndx, end_vertex_ndx);
      this.addEdgeBetweenNodes(start_node_ndx, end_node_ndx, min_gain, max_gain);
    end % End of function
  end % End of protected methods.

  methods(Access = {?matlab.unittest.TestCase, ?FlowTransitionGainDigraph}) % Private block, accessible by unit tests.
    
    % ╭─────────────────────────────────────────────────────────────╮
    % │  ╭───────────────────────────────────────────────────────╮  │
    % │  │             Operations Using Node Indices             │  │
    % │  ╰───────────────────────────────────────────────────────╯  │
    % ╰─────────────────────────────────────────────────────────────╯
    function [node_row] = getNodeRow(this, node_ndx)
      arguments(Input)
        this;
        node_ndx; 
      end % End of Input arguments block
      node_row = this.gains_digraph.Nodes(node_ndx, :);
    end % End of function

    function edge_rows = getEdgesBetweenNodes(this, start_node_ndx, end_node_ndx)
      arguments(Input)
        this;
        start_node_ndx {pwintz.validators.mustBeIndexScalar};
        end_node_ndx   {pwintz.validators.mustBeIndexScalar};
      end % End of Input arguments block
      edge_ndxs = this.getEdgeIndicesBetweenNodes(start_node_ndx, end_node_ndx);
      edge_rows = this.gains_digraph.Edges(edge_ndxs, :);
    end % End of function

    function edge_ndxs = getEdgeIndicesBetweenNodes(this, start_node_ndx, end_node_ndx)
      arguments(Input)
        this;
        start_node_ndx {pwintz.validators.mustBeIndexScalar};
        end_node_ndx   {pwintz.validators.mustBeIndexScalar};
      end % End of Input arguments block
      
      edge_to_find = [start_node_ndx, end_node_ndx];
      edge_ndxs = pwintz.arrays.findRowIn(edge_to_find, this.gains_digraph.Edges.EndNodes);
    end % End of function

    function has_edge = hasEdgeBetweenNodes(this, start_node_ndx, end_node_ndx)
      has_edge = ~isempty(this.getEdgeIndicesBetweenNodes(start_node_ndx, end_node_ndx));
    end % End of function

    function this = addEdgeBetweenNodes(this, start_node_ndx, end_node_ndx, min_gain, max_gain)
      assert(start_node_ndx <= this.numNodes());
      assert(start_node_ndx >= 1);
      assert(end_node_ndx >= 1);
      assert(min_gain >= 0);
      assert(max_gain >= 0);
      assert(min_gain <= max_gain);

      start_name = this.gains_digraph.Nodes.Names(start_node_ndx);
      end_name   = this.gains_digraph.Nodes.Names(end_node_ndx);

      edge_table = struct2table(struct(...
        "EndNodes", [start_node_ndx, end_node_ndx], ...
        "Label",   pwintz.strings.format("(%s) -> (%s)", start_name, end_name), ...
        "MaxGain",  max_gain, ...
        "MinGain",  min_gain ...
      ));
      this.gains_digraph = this.gains_digraph.addedge(edge_table);
    end % End of function

    function edge_ndxs = getEdgeIndicesFromVertexToCone(this, start_vertex_ndx, end_cone_ndx)
      start_node_ndx = this.nodeIndexOfModeAndVertex(start_vertex_ndx);
      end_node_ndx = this.nodeIndexOfModeAndCone(end_cone_ndx);
      edge_ndxs = this.getEdgeIndicesBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function


    % ╭───────────────────────────────╮
    % │             Paths             │
    % ╰───────────────────────────────╯
    function [paths_node_ndxs, paths_edge_ndxs, paths_edge_rows] = getPathsBetweenNodes(this, start_node_ndx, end_node_ndx)
      % We use MinPathLength=1 to prevent allpaths from producing an empty path when start_vertex_ndx == end_vertex_ndx.
      [paths_node_ndxs, paths_edge_ndxs] = this.gains_digraph.allpaths(start_node_ndx, end_node_ndx, 'MinPathLength', 1);
      paths_edge_rows = cellfun(@(edge_ndxs) this.gains_digraph.Edges(edge_ndxs, :), paths_edge_ndxs, 'UniformOutput', false);
    end % End of function


    % ╭──────────────────────────────────────╮
    % │             Node Indices             │
    % ╰──────────────────────────────────────╯

    function node_ndx = nodeIndexOfModeAndVertex(this, mode_ndx, vertex_ndx)
      arguments(Input)
        this;
        mode_ndx   (1, 1) {pwintz.validators.mustBeIndexScalar};
        vertex_ndx (:, 1) {pwintz.validators.mustBeIndexVector};
      end % End of Input arguments block

      pwintz.assertions.assertAllAreMembers(mode_ndx,   this.mode_indices, arrayName="mode_indices");
      pwintz.assertions.assertAllAreMembers(vertex_ndx, this.conical_partitions{mode_ndx}.vertex_indices, arrayName="vertex_indices");

      vertex_node_ndxs_in_mode = this.node_ndxs_of_mode_and_vertices{mode_ndx};
      node_ndx = vertex_node_ndxs_in_mode(vertex_ndx);
      
      node_row = this.gains_digraph.Nodes(node_ndx, :);
      pwintz.assertions.assertEqual(node_row.ConeIndex, 0);
    end
    
    function node_ndx = nodeIndexOfModeAndCone(this, mode_ndx, cone_ndx)
      arguments(Input)
        this;
        mode_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        cone_ndx (:, 1) {pwintz.validators.mustBeIndexVector};
      end % End of Input arguments block

      pwintz.assertions.assertAllAreMembers(mode_ndx, this.mode_indices);
      pwintz.assertions.assertAllAreMembers(cone_ndx, this.conical_partitions{mode_ndx}.cone_indices, arrayName="cone_indices");

      cone_node_ndxs_in_mode = this.node_ndxs_of_mode_and_cones{mode_ndx};
      node_ndx = cone_node_ndxs_in_mode(cone_ndx);

      node_row = this.gains_digraph.Nodes(node_ndx, :);
      pwintz.assertions.assertEqual(node_row.VertexIndex, 0);
    end

    function node_ndx = nodeIndexOfVertex(this, vertex_ndx)
      % assert(isscalar(vertex_ndx), "nodeIndexOfVertex(): vertex_ndx must be a scalar.");
      % assert(vertex_ndx >= 0, "Invalid vertex index '%s'. Must be in {1, 2, ..., %d}", mat2str(vertex_ndx), this.numVertexNodes);
      % assert(vertex_ndx <= this.numVertexNodes, "Invalid vertex index '%s'. Must be in {1, 2, ..., %d}", mat2str(vertex_ndx), this.numVertexNodes);
      % node_ndx = this.node_ndxs_of_mode_and_vertices(vertex_ndx);
      error("Deprecated. Use nodeIndexOfModeAndVertex");
    end % End of function
    
    function node_ndx = nodeIndexOfCone(this, cone_ndx)
      % assert(isscalar(cone_ndx), "nodeIndexOfCone(): cone_ndx must be a scalar.");
      pwintz.assertions.assertAllGreaterThanOrEqual(cone_ndx, 0,           "Invalid cone index '%s'. Must be in {1, 2, ..., %d}", mat2str(cone_ndx), this.numConeNodes);
      pwintz.assertions.assertAllLessThanOrEqual(cone_ndx, this.numConeNodes, "Invalid cone index '%s'. Must be in {1, 2, ..., %d}", mat2str(cone_ndx), this.numConeNodes);
      node_ndx = this.node_ndxs_of_mode_and_cones(cone_ndx);
      error("Deprecated. Use nodeIndexOfModeAndCone");
    end % End of function

    function vertex_ndx = vertexIndexOfNode(this, node_ndx_of_vertex)
      arguments(Input)
        this;
        node_ndx_of_vertex (1, 1) {mustBeScalarOrEmpty,mustBeNonempty,mustBeInteger,mustBePositive};
      end % End of Input arguments block
      assert(node_ndx_of_vertex <= this.numVertexNodes(), "vertexIndexOfNode(): Invalid vertex node index '%s'. Must be in <=  number of vertex nodes.", mat2str(node_ndx_of_vertex), this.numVertexNodes());
      vertex_ndx = node_ndx_of_vertex;
    end % End of function

    function cone_ndx = coneIndexOfNodeIndex(this, node_ndx_of_cone)
      arguments(Input)
        this;
        node_ndx_of_cone (1, 1) {mustBeScalarOrEmpty,mustBeNonempty,mustBeInteger,mustBePositive};
      end % End of Input arguments block
      assert(node_ndx_of_cone >= this.numVertexNodes(), "coneIndexOfNodeIndex(): Invalid vertex node index '%s'. Must be in > number of vertices.", mat2str(node_ndx_of_cone), this.numVertexNodes());
      assert(node_ndx_of_cone <= this.numNodes(), "coneIndexOfNodeIndex(): Invalid vertex node index '%s'. Must be in <= number of nodes.", mat2str(node_ndx_of_cone), this.numNodes());
      cone_ndx = node_ndx_of_cone - this.numVertexNodes();
    end % End of function
    
  end % End private methods block
end % End of class

function label_array = coneIndicesToTexLabels(node_ndx, cone_ndxs)
  arguments(Input)
    node_ndx (1, 1) double {mustBeInteger,mustBePositive};
    cone_ndxs (1, :) double {mustBeInteger,mustBePositive};
  end % End of Input arguments block
  arguments(Output)
    label_array (:, 1) string;% cell;
  end % End of Output arguments block
  % label_array = cellfun(@(ndx) char(sprintf("$C_{%d}", ndx)), cone_ndxs, 'UniformOutput', false);
  label_array = arrayfun(@(ndx) sprintf("$C^{%d}_{%d}$", node_ndx, ndx), cone_ndxs, 'UniformOutput', false);
end

function label_array = vertexIndicesToTexLabels(node_ndx, vertex_ndxs)
  arguments(Input)
    node_ndx (1, 1) double {mustBeInteger,mustBePositive};
    vertex_ndxs (1, :) double {mustBeInteger,mustBePositive};
  end % End of Input arguments block
  arguments(Output)
    label_array (:, 1) string;
  end % End of Output arguments block
  label_array = arrayfun(@(ndx) sprintf("$v^{%d}_{%d}$", node_ndx, ndx), vertex_ndxs, 'UniformOutput', false);
end % end function