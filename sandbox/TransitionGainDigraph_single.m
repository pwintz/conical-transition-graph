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
    node_ndxs_of_vertices (:,1); % Row vector of natural numbers.
    node_ndxs_of_cones    (:,1); % Row vector of natural numbers.
    conical_partition     (1,1); % ConicalPartition;
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
    function union_tgd = union(tgd_left, tgd_right)
      pwintz.assertions.assertEqual(tgd_left.gains_digraph.Nodes, tgd_right.gains_digraph.Nodes);
      
      union_tgd = TransitionGainDigraph(tgd_left.conical_partition);
      union_edges = [
        tgd_left.gains_digraph.Edges;
        tgd_right.gains_digraph.Edges;
      ];
      % union_tgd.gains_digraph.Edges = union_edges;
      union_tgd.gains_digraph = union_tgd.gains_digraph.addedge(union_edges);
    end % End of function
  end % End static methods block

  methods
    
    % Constructor
    function this = TransitionGainDigraph(conical_partition)
      Names = ["Vertex " + conical_partition.vertex_indices'; "Cone " + conical_partition.cone_indices'];

      VertexIndex = [conical_partition.vertex_indices'; 0 * conical_partition.cone_indices'];
      ConeIndex   = [0*conical_partition.vertex_indices'; conical_partition.cone_indices'];

      this.node_ndxs_of_vertices = find(VertexIndex);
      this.node_ndxs_of_cones    = find(ConeIndex);

      VertexTeXLabels = vertexIndicesToTexLabels(conical_partition.vertex_indices);
      ConeTeXLabels   = coneIndicesToTexLabels(conical_partition.cone_indices);

      Objects = [num2cell(conical_partition.vertices, 1), num2cell(conical_partition.cones)]';
      TeXLabel = [VertexTeXLabels; ConeTeXLabels];

      cone_centers = arrayfun(@(cone) centroid(cone.rays')', conical_partition.cones, "UniformOutput", false);
      % cone_centers = cell2mat(cellfun(@(convex_poly) centroid(convex_poly.vertices')', conical_partition.cones, "UniformOutput", false));
      vert_cell_array = num2cell(conical_partition.vertices, 1);
      PlotPosition = [vert_cell_array, cone_centers];

      % ╭────────────────────────────────────────────────╮
      % │             Construct gain_digraph             │
      % ╰────────────────────────────────────────────────╯
      if exist("gain_digraph", "var")
        this.gains_digraph = gains_digraph;
      else
        PlotPosition = [PlotPosition{:}]';        
        node_table = pwintz.tables.makeTable(...
          "Names", Names, ...
          "Objects", Objects, ...
          "PlotPosition", PlotPosition, ...
          "VertexIndex", VertexIndex, ...
          "ConeIndex", ConeIndex, ...
          "TeXLabel", TeXLabel...
        );

        % node_table = table(...
        %   Names,          ...
        %   Objects,        ...
        %   PlotPosition,   ...
        %   VertexIndex,    ...
        %   ConeIndex,      ...
        %   TeXLabel        ...
        % )
        % node_table

        edge_table = struct2table(struct(...
          "EndNodes", double.empty(0, 2), ...
          ..."Weight",   double.empty(0, 1), ...
          "MaxGain",  double.empty(0, 1), ...
          "MinGain",  double.empty(0, 1)  ...
        ));
        this.gains_digraph = digraph(edge_table, node_table);
      end % End of if block.
      this.conical_partition = conical_partition;
    end
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
    function n = numVertexNodes(this)
      n = numel(this.node_ndxs_of_vertices);
      % n = numel(this.gains_digraph.Edges);
    end % End of function
    function n = numConeNodes(this)
      n = numel(this.node_ndxs_of_cones);
    end % End of function

    % ╭────────────────────────────────────────────────╮
    % │             Get Items By Reference             │
    % ╰────────────────────────────────────────────────╯
    function vertex_row = getVertexRow(this, vertex_ndx)
      node_ndx = this.nodeIndexOfVertex(vertex_ndx);
      vertex_row = this.gains_digraph.Nodes(node_ndx, :);
    end % End of function

    function cone_row = getConeRow(this, cone_ndx)
      node_ndx = this.nodeIndexOfCone(cone_ndx);
      cone_row = this.gains_digraph.Nodes(node_ndx, :);
    end % End of function

    function edge_row = getEdgeRow(this, edge_ndx)
      edge_row = this.gains_digraph.Edges(edge_ndx, :);
    end % End of function

    % ╭──────────────────────────────────────────╮
    % │             Has Edge Methods             │
    % ╰──────────────────────────────────────────╯
    function has_edge = hasEdgeFromVertexToCone(this, start_vertex_ndx, end_cone_ndx)
      start_node_ndx = this.nodeIndexOfVertex(start_vertex_ndx);
      end_node_ndx = this.nodeIndexOfCone(end_cone_ndx);
      has_edge = this.hasEdgeBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function

    function has_edge = hasEdgeFromVertexToVertex(this, start_vertex_ndx, end_vertex_ndx)
      start_node_ndx = this.nodeIndexOfVertex(start_vertex_ndx);
      end_node_ndx = this.nodeIndexOfVertex(end_vertex_ndx);
      has_edge = this.hasEdgeBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function

    function has_edge = hasEdgeFromConeToCone(this, start_cone_ndx, end_cone_ndx)
      start_node_ndx = this.nodeIndexOfCone(start_cone_ndx);
      end_node_ndx = this.nodeIndexOfCone(end_cone_ndx);
      has_edge = this.hasEdgeBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function

    function has_edge = hasEdgeFromConeToVertex(this, start_cone_ndx, end_vertex_ndx)
      start_node_ndx = this.nodeIndexOfCone(start_cone_ndx);
      end_node_ndx = this.nodeIndexOfVertex(end_vertex_ndx);
      has_edge = this.hasEdgeBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function

    % ╭─────────────────────────────────────────╮
    % │             Get Graph Edges             │
    % ╰─────────────────────────────────────────╯
    function edge_rows = getEdgesFromVertexToCone(this, start_vertex_ndx, end_cone_ndx)
      arguments(Input)
        this;
        start_vertex_ndx (:, 1) {pwintz.validators.mustBeIndexVector};
        end_cone_ndx     (:, 1) {pwintz.validators.mustBeIndexVector};
      end % End of Input arguments block
      
      assert(all(start_vertex_ndx <= size(this.node_ndxs_of_vertices, 1)), "start_vertex_ndx=%s must be less than %d", mat2str(start_vertex_ndx), size(this.node_ndxs_of_vertices, 1));
      assert(all(end_cone_ndx <= size(this.node_ndxs_of_cones, 1)), "end_cone_ndx=%s must be less than %d", mat2str(end_cone_ndx), size(this.node_ndxs_of_cones, 1));
      start_node_ndx = this.node_ndxs_of_vertices(start_vertex_ndx);
      end_node_ndx   = this.node_ndxs_of_cones(end_cone_ndx);
      edge_rows      = this.getEdgesBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function
    
    function edge_rows = getEdgesFromConeToVertex(this, start_cone_ndx, end_vertex_ndx)
      start_node_ndx = this.node_ndxs_of_cones(start_cone_ndx);
      end_node_ndx   = this.node_ndxs_of_vertices(end_vertex_ndx);
      edge_rows      = this.getEdgesBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function
    
    function edge_rows = getEdgesFromConeToCone(this, start_cone_ndx, end_cone_ndx)    
      start_node_ndx = this.node_ndxs_of_cones(start_cone_ndx);
      end_node_ndx   = this.node_ndxs_of_cones(end_cone_ndx);  
      edge_rows      = this.getEdgesBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function
    
    function edge_rows = getEdgesFromVertexToVertex(this, start_vertex_ndx, end_vertex_ndx)
      start_node_ndx = this.node_ndxs_of_vertices(start_vertex_ndx);
      end_node_ndx   = this.node_ndxs_of_vertices(end_vertex_ndx);
      edge_rows      = this.getEdgesBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function

    % ╭───────────────────────────────────────╮
    % │         All Edges, by Type            │
    % ╰───────────────────────────────────────╯
    function edge_rows = getEdgesFromVerticesToVertices(this)
      all_edges = this.gains_digraph.Edges;
      edge_start_node_ndxs = all_edges.EndNodes(:, 1);
      edge_end_node_ndxs   = all_edges.EndNodes(:, 2);
      
      does_edge_start_at_vertex = ismember(edge_start_node_ndxs, this.node_ndxs_of_vertices);
      does_edge_end_at_vertex   = ismember(edge_end_node_ndxs,   this.node_ndxs_of_vertices);

      edge_rows = all_edges(does_edge_start_at_vertex & does_edge_end_at_vertex, :);
      
    end % End of function

    function edge_rows = getEdgesFromVerticesToCones(this)
      all_edges = this.gains_digraph.Edges;
      edge_start_node_ndxs = all_edges.EndNodes(:, 1);
      edge_end_node_ndxs   = all_edges.EndNodes(:, 2);
      
      does_edge_start_at_vertex = ismember(edge_start_node_ndxs, this.node_ndxs_of_vertices);
      does_edge_end_at_cone   = ismember(edge_end_node_ndxs,   this.node_ndxs_of_cones);
      
      edge_rows = all_edges(does_edge_start_at_vertex & does_edge_end_at_cone, :);
    end % End of function
    
    function edge_rows = getEdgesFromConesToVertices(this)
      all_edges = this.gains_digraph.Edges;
      edge_start_node_ndxs = all_edges.EndNodes(:, 1);
      edge_end_node_ndxs   = all_edges.EndNodes(:, 2);
      
      does_edge_start_at_cone = ismember(edge_start_node_ndxs, this.node_ndxs_of_cones);
      does_edge_end_at_vertex   = ismember(edge_end_node_ndxs,   this.node_ndxs_of_vertices);
      
      edge_rows = all_edges(does_edge_start_at_cone & does_edge_end_at_vertex, :);
    end % End of function


    function edge_rows = getEdgesFromConesToCones(this)
      all_edges = this.gains_digraph.Edges;
      edge_start_node_ndxs = all_edges.EndNodes(:, 1);
      edge_end_node_ndxs   = all_edges.EndNodes(:, 2);
      
      does_edge_start_at_cone = ismember(edge_start_node_ndxs, this.node_ndxs_of_cones);
      does_edge_end_at_cone   = ismember(edge_end_node_ndxs,   this.node_ndxs_of_cones);
      
      edge_rows = all_edges(does_edge_start_at_cone & does_edge_end_at_cone, :);
      assert(isempty(edge_rows));
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
      start_node_ndx = this.node_ndxs_of_vertices(start_vertex_ndx);
      end_node_ndx   = this.node_ndxs_of_vertices(end_vertex_ndx);
      [paths_node_ndxs, paths_edge_ndxs, paths_edge_rows] = this.getPathsBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function

    function [paths_node_ndxs, paths_edge_ndxs, paths_edge_rows] = getPathsBetweenCones(this, start_cone_ndx, end_cone_ndx)
      start_node_ndx = this.node_ndxs_of_cones(start_cone_ndx);
      end_node_ndx   = this.node_ndxs_of_cones(end_cone_ndx);
      [paths_node_ndxs, paths_edge_ndxs, paths_edge_rows] = this.getPathsBetweenNodes(start_node_ndx, end_node_ndx);
    end % End of function

    function [cycles_nodes, cycles_edges] = getCycles(this)
      [cycles_nodes, cycles_edges] = this.gains_digraph.allcycles();
    end % End of function

    function [cycles_nodes, cycles_edges] = getVertexCycles(this)
      vertex_graph = subgraph(this.gains_digraph, this.node_ndxs_of_vertices)  ;
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
    function plot(this)
      graph = this.gains_digraph;
      positions = graph.Nodes.PlotPosition;
      p = plot(graph, ...
        "XData", positions(:, 1), ...
        "YData", positions(:, 2), ...
        "NodeColor", "black", ...
        "EdgeAlpha", 0.5, ...
        "MarkerSize", 8, ...
        "LineWidth", 1, ...
        ...'NodeLabel', graph.Nodes.TeXLabel, ...
        ...'NodeLabelMode', 'manual', ...
        ...'EdgeLabel', graph.Edges.WeightIntervalStr, ...
        "NodeFontSize", 12, ...
        ... "EdgeFontSize", 10 * graph.Edges.Length / max(graph.Edges.Length), ...
        "Interpreter", "latex"...
      );
      
      % !! Using the TeX labels does not work for unknown reasons. The result is no label at all, so for now we will just use the default "node number" labels.
      ...p.NodeLabel = graph.Nodes.TeXLabel;

      if this.numEdges() > 0
        geq_1_color = [1 0 0];
        lt_1_color  = [0 0 1];
        edge_color  = (graph.Edges.MaxGain >= 1) * geq_1_color + (graph.Edges.MaxGain < 1) * lt_1_color;
        p.EdgeColor = edge_color;
        
        % For the plotted weights, we want to illustrate that the importance aspect of an edges weight is its distance above or below 1. We use abs(log(weight)) to normalize. This way, weights equal to 1/2 and 2 are plotted with the same width.
        weight = abs(log(graph.Edges.MaxGain));
        weight = 6 * (weight + 1) / (max(weight) + 1) + 1;
        % p.LineWidth = weight;
        p.ArrowSize = 8*sqrt(p.LineWidth);
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
%       ismember(this.gains_digraph.Edges.EndNodes, this.node_ndxs_of_cones)
% 
%       this.gains_digraph.rmedge()
%       vertex_subgraph = this.gains_digraph.subgraph(this.node_ndxs_of_vertices);
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
%       % vertex_subgraph = this.gains_digraph.subgraph(this.node_ndxs_of_vertices);
%     end % End of function
% 
%     function vertex_subgraph = getConeSubgraph(this)
%       arguments(Output)
%         vertex_subgraph (1, 1) digraph;
%       end % End of Output arguments block
%       
%       vertex_subgraph = this.gains_digraph.subgraph(this.node_ndxs_of_cones);
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
      string_rep = sprintf("%s with %d nodes (%d vertex nodes and %d cone nodes) and %d edges.", class(this), this.numNodes(), this.numVertexNodes(), this.numConeNodes(), this.numEdges());
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


  methods(Access=protected) % Define protected methods.
    % ╭─────────────────────────────────────────╮
    % │  ╭───────────────────────────────────╮  │
    % │  │             Add Edges             │  │
    % │  ╰───────────────────────────────────╯  │
    % ╰─────────────────────────────────────────╯
    function this = addEdgeFromVertexToCone(this, start_vertex_ndx, end_cone_ndx, min_gain, max_gain)
      arguments(Input)
        this;
        start_vertex_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        end_cone_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        min_gain (1, 1) double {mustBeNonnegative};
        max_gain (1, 1) double {mustBeNonnegative};
      end % End of Input arguments block
      
      start_node_ndx = this.node_ndxs_of_vertices(start_vertex_ndx);
      end_node_ndx   = this.node_ndxs_of_cones(end_cone_ndx);
      this.addEdgeBetweenNodes(start_node_ndx, end_node_ndx, min_gain, max_gain);
    end % End of function
  
    function this = addEdgeFromConeToVertex(this, start_cone_ndx, end_vertex_ndx, min_gain, max_gain)
      arguments(Input)
        this;
        start_cone_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        end_vertex_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        min_gain (1, 1) double {mustBeNonnegative};
        max_gain (1, 1) double {mustBeNonnegative};
      end % End of Input arguments block
      
      start_node_ndx = this.node_ndxs_of_cones(start_cone_ndx);
      end_node_ndx   = this.node_ndxs_of_vertices(end_vertex_ndx);
      this.addEdgeBetweenNodes(start_node_ndx, end_node_ndx, min_gain, max_gain);
    end % End of function
  
    function this = addEdgeFromConeToCone(this, start_cone_ndx, end_cone_ndx, min_gain, max_gain)    
      arguments(Input)
        this;
        start_cone_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        end_cone_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        min_gain (1, 1) double {mustBeNonnegative};
        max_gain (1, 1) double {mustBeNonnegative};
      end % End of Input arguments block
      start_node_ndx = this.node_ndxs_of_cones(start_cone_ndx);
      end_node_ndx   = this.node_ndxs_of_cones(end_cone_ndx);  
      this.addEdgeBetweenNodes(start_node_ndx, end_node_ndx, min_gain, max_gain);
    end % End of function
  
    function this = addEdgeFromVertexToVertex(this, start_vertex_ndx, end_vertex_ndx, min_gain, max_gain)
      arguments(Input)
        this;
        start_vertex_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        end_vertex_ndx (1, 1) {pwintz.validators.mustBeIndexScalar};
        min_gain (1, 1) double {mustBeNonnegative};
        max_gain (1, 1) double {mustBeNonnegative};
      end % End of Input arguments block
      start_node_ndx = this.node_ndxs_of_vertices(start_vertex_ndx);
      end_node_ndx   = this.node_ndxs_of_vertices(end_vertex_ndx);
      this.addEdgeBetweenNodes(start_node_ndx, end_node_ndx, min_gain, max_gain);
    end % End of function
  end % End of protected methods.

  methods(Access = {?matlab.unittest.TestCase}) % Private block, accessible by unit tests.
    
    % ╭─────────────────────────────────────────────────────────────╮
    % │  ╭───────────────────────────────────────────────────────╮  │
    % │  │             Operations Using Node Indices             │  │
    % │  ╰───────────────────────────────────────────────────────╯  │
    % ╰─────────────────────────────────────────────────────────────╯
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
      assert(start_node_ndx <= this.numNodes());
      assert(start_node_ndx >= 1);
      assert(end_node_ndx >= 1);
      assert(min_gain >= 0);
      assert(max_gain >= 0);
      assert(max_gain >= min_gain);
      weight = max_gain;
      edge_table = struct2table(struct(...
        "EndNodes", [start_node_ndx, end_node_ndx], ...
        ..."Weight",   weight, ...
        "MaxGain",  max_gain, ...
        "MinGain",  min_gain ...
      ));
      this.gains_digraph = this.gains_digraph.addedge(edge_table);
    end % End of function

    function edge_ndxs = getEdgeIndicesFromVertexToCone(this, start_vertex_ndx, end_cone_ndx)
      start_node_ndx = this.nodeIndexOfVertex(start_vertex_ndx);
      end_node_ndx = this.nodeIndexOfCone(end_cone_ndx);
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
    function node_ndx = nodeIndexOfVertex(this, vertex_ndx)
      assert(isscalar(vertex_ndx), "nodeIndexOfVertex(): vertex_ndx must be a scalar.");
      assert(vertex_ndx >= 0, "Invalid vertex index '%s'. Must be in {1, 2, ..., %d}", mat2str(vertex_ndx), this.numVertexNodes);
      assert(vertex_ndx <= this.numVertexNodes, "Invalid vertex index '%s'. Must be in {1, 2, ..., %d}", mat2str(vertex_ndx), this.numVertexNodes);
      node_ndx = this.node_ndxs_of_vertices(vertex_ndx);
    end % End of function
    
    function node_ndx = nodeIndexOfCone(this, cone_ndx)
      % assert(isscalar(cone_ndx), "nodeIndexOfCone(): cone_ndx must be a scalar.");
      pwintz.assertions.assertAllGreaterThanOrEqual(cone_ndx, 0,           "Invalid cone index '%s'. Must be in {1, 2, ..., %d}", mat2str(cone_ndx), this.numConeNodes);
      pwintz.assertions.assertAllLessThanOrEqual(cone_ndx, this.numConeNodes, "Invalid cone index '%s'. Must be in {1, 2, ..., %d}", mat2str(cone_ndx), this.numConeNodes);
      node_ndx = this.node_ndxs_of_cones(cone_ndx);
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

function label_array = coneIndicesToTexLabels(cone_ndxs)
  arguments(Input)
    cone_ndxs (1, :) double {mustBeInteger,mustBePositive};
  end % End of Input arguments block
  arguments(Output)
    label_array (:, 1) string;% cell;
  end % End of Output arguments block
  % label_array = cellfun(@(ndx) char(sprintf("$C_{%d}", ndx)), cone_ndxs, 'UniformOutput', false);
  label_array = arrayfun(@(ndx) sprintf("$C_{%d}", ndx), cone_ndxs, 'UniformOutput', false);
end

function label_array = vertexIndicesToTexLabels(vertex_ndxs)
  arguments(Input)
    vertex_ndxs (1, :) double {mustBeInteger,mustBePositive};
  end % End of Input arguments block
  arguments(Output)
    label_array (:, 1) string;
  end % End of Output arguments block
  label_array = arrayfun(@(ndx) sprintf("$v_{%d}", ndx), vertex_ndxs, 'UniformOutput', false);
end % end function