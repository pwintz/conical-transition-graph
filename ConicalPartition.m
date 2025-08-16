classdef ConicalPartition < handle
  %%% ConicalPartition - Create a class for defining a partition of R^n into cones. 
  %%% The cones are defined usings a "watermelon slicing" method using equal angles for all cones.
  %%% To keep things simple, we are starting with  

  %%%
  
  properties(SetAccess = immutable)
    % Define instance constants.
    A % Matrix of "\dot x = A x".
    dimension  % The dimension of the system considered
    derivative_vertices
    state_vertices
  end
  properties(SetAccess = immutable, GetAccess=private)
    nSlices % Number of divisions per dimension  
    derivativeSpaceAngle % The angle swept by each cone in derivative space along the rotation in each dimension.
    derivativeSpaceSliceAngles % The angle of each slice in derivative space.
    stateVertexAngles
    vertexMultiIndexSize
  end

  properties(Dependent)
    % Define dependent variables.
    n_state_vertices
    state_vertex_indices
  end
  methods
    function value = get.n_state_vertices(this)
      value = size(this.state_vertices, 2);
    end
    function value = get.state_vertex_indices(this)
      value = 1:this.n_state_vertices;
    end
  end

  methods
    % Constructor
    function this = ConicalPartition(varargin)
      % Parse inputs from varargin argument.
      p = inputParser;
      addRequired(p, "A", @(x) isnumeric(x) && ismatrix(x)); % Positional argument.
      addParameter(p, "dimension", 2, @(x) x == 2); % Name-value parameter.
      addParameter(p, "nSlices", 8, @(x) isreal(x) && rem(x,1)==0 && x > 1); % Name-value parameter.;
      
      parse(p, varargin{:});
      args = p.Results;
      this.A = args.A;
      this.dimension    = args.dimension;
      this.nSlices      = args.nSlices;
      this.vertexMultiIndexSize = 1; % For 2D Case. Must be assigned before the first use of get...Vertex()

      % ⋘──── Compute Property Values Derived from Given Values ────⋙
      angles = linspace(0, 2*pi, this.nSlices + 1);
      this.derivativeSpaceSliceAngles = angles(1:end-1);
      this.derivative_vertices = nan([2, this.nSlices]);
      this.state_vertices = nan([2, this.nSlices]);
      % this.derivativeVertexAngles = nan([1, this.nSlices]);
      this.stateVertexAngles = nan([1, this.nSlices]);
      for i = 1:this.nSlices
        this.derivative_vertices(:,i) = this.getDerivativeSpaceVertex(i);
        this.state_vertices(:,i) = this.getStateSpaceVertex(i);
        this.stateVertexAngles(i) = atan2(this.state_vertices(2,i), this.state_vertices(1,i));
      end
      % ⋘────────── Check Properties   ──────────⋙
      assert(size(this.A,1) == size(this.A,2), "The matrix A must be square")
      assert(cond(this.A) < 1e9, "The matrix A must be well-conditioned but had a condition number of %8.2g", cond(this.A))
    end

    % ╭──────────────────────────────────╮
    % │             Indexing             │
    % ╰──────────────────────────────────╯
    function multiIndexArray = getVertexMultiIndices(this)
      % vertexMultiIndices - Get an wide array where each column is the multi-index of a one vertex in this conical partition.
      % The result can be used like 
      % *
      % *  for i = conical_partition.getVertexMultiIndices()
      % *    v_i = conical_partition.getStateSpaceVertex(i)
      % *  end
      % *

      % For the 2D case, each multiIndexing only has one entry.
      multiIndexArray = 1:this.nSlices;
    end
    
    function v_neighbors_ndxs = getNeighborVertexIndices(this, multiIndex)
      % getNeighborVertexIndices - Get an wide array containing where each column is the multi-index of a vertex that is adjacent to the given multi-index.
      % The result can be used like 
      % *                                                
      % * for i = conical_partition.getVertexMultiIndices() 
      % *   v_i = conical_partition.getStateSpaceVertex(i)         
      % *   for j = conical_partition.getNeighborVertexIndices(i) 
      % *     v_j = conical_partition.getStateSpaceVertex(i)  
      % *   end                                                   
      % * end                                                   
      % *                                                

      neighborKernel = [-1, +1]; % For 2D case

      v_neighbors_ndxs = this.normalizeIndex(multiIndex + neighborKernel);
    end

     function adjacent_regions_ndxs = getAdjacentRegionsIndices(this, v_ndx)
      % getAdjacentRegionsIndices - Get a list of multi-indices for regions with the given vertex (defined by v_ndx) in their boundary.
      % The result is wide array where each column is the multi-index of a region.
      neighborRegionsNdxOffsets = [-1, 0]; % For 2D case

      adjacent_regions_ndxs = this.normalizeIndex(v_ndx + neighborRegionsNdxOffsets);
    end

    function conjoining_regions_ndxs = getConjoiningRegionsIndices(this, v1_ndx, v2_ndx)
      % getConjoiningRegionsIndices - Get a list of regions that are adjacent to two given vertices.
      % The result is wide array where each column is the multi-index of a region.
      v1_adjacent_regions_ndxs = this.getAdjacentRegionsIndices(v1_ndx);
      v2_adjacent_regions_ndxs = this.getAdjacentRegionsIndices(v2_ndx);

      
      conjoining_regions_ndxs = intersect(v1_adjacent_regions_ndxs, v2_adjacent_regions_ndxs);
    end

    function v = getStateSpaceVertex(this, multiIndex)
      % getVertex - Given the multi-index of the vertex, passed as the arguments, returns the vertex as a column vector.
      vertex_index = this.normalizeIndex(multiIndex);
      v = this.getDerivativeSpaceVertex(vertex_index);
      v = inv(this.A)*v;
      v = v / norm(v);
    end

    function v = getDerivativeSpaceVertex(this, multiIndex)
      % getVertex - Given the multi-index of the vertex, passed as the arguments, returns the vertex as a column vector.
      vertex_index = this.normalizeIndex(multiIndex);
      vertex_angle = this.derivativeSpaceSliceAngles(vertex_index);
      v = [cos(vertex_angle); sin(vertex_angle)];

    end

    function theta = getStateSpaceVertexAngle(this, multiIndex)
      % getVertex - Given the multi-index of the vertex, passed as the arguments, returns the vertex as a column vector.
      vertex_index = this.normalizeIndex(multiIndex);
      theta = this.stateVertexAngles(vertex_index);
    end
    
    function C = getStateCone(this, multiIndex)
      % Given an index, get the state space cone for a given multi-index, given as a polyhedron (the cone truncated to the unit sphere), which is encoded as a row of vertices that defines the polyhedron.
      vertex_index = this.normalizeIndex(multiIndex);
      v1 = this.getStateSpaceVertex(vertex_index);
      v2 = this.getStateSpaceVertex(vertex_index + 1);
      origin = [0; 0];
      C = ConvexPolyhedron.fromConvexHull([v1, origin, v2]);
    end

    function D = getDerivativeCone(this, multiIndex)
      % Given an index, get the state space cone for a given multi-index, given as a polyhedron (the cone truncated to the unit sphere), which is encoded as a row of vertices that defines the polyhedron.
      vertex_index = this.normalizeIndex(multiIndex);
      v1 = this.getDerivativeSpaceVertex(vertex_index);
      v2 = this.getDerivativeSpaceVertex(vertex_index+1);
      origin = [0; 0];
      D = ConvexPolyhedron.fromConvexHull([v1, origin, v2]);
    end
    
    % ╭────────────────────────────────────────╮
    % │             Reachable Sets             │
    % ╰────────────────────────────────────────╯
    function R = computeReachableSet(this, regionMultiIndex, P0)
      % Compute the region that is reachable from P0 when flowing according to \dot x \in AD, where D is the derivative cone for the given multiIndex.
      D = conical_partition.getDerivativeCone(regionMultiIndex);
      % The "derivative cone" is represented using vertices on the unit sphere, so we scale it up so that we get a set the sufficiently represents "P0 + AD", at least up to a large radius.
      R = addDerivativeCone(P, 1000*D);
      C = conical_partition.getStateCone(i);
      R_in_C = intersection(R, 1000*C);
    end
  end

  methods(Access=private) % Define private class methods.

    function out = normalizeIndex(this, multiIndex)
      % normalizeIndex - Given a multi-index, normalize its entries by applying the modulo operation so that each entry is in {1, 2, ..., nSlices}.

      % Check that the multi-index matches the dimension of the system.
      assert( ...
        size(multiIndex, 1) == this.vertexMultiIndexSize, ...
        "The size(multiIndex, 1) = %d must be vertexMultiIndexSize = %d.", size(multiIndex, 1), this.vertexMultiIndexSize ...
      )

      % Compute the index of the vertex from the given argument, modulo the number of slices. Since MATLAB uses 1-based indexing, we need to subtract 1 before computing the modulo and then add it afterward.
      out = mod(multiIndex - 1, this.nSlices) + 1;
    end
  end
end

function result = addDerivativeCone(P, D)
  % Given a convex polygon P, add the derivative cone D, which consists of three points [v1, 0, v2]. 
  % The addition is in the sense of the Minkowski sum.
  v1 = D(:, 1);
  v2 = D(:, 3);
  result = [P, P + v1, P + v2];
  dTriangulated = delaunayTriangulation(result');
  convexHullIndices = convexHull(dTriangulated);
  result = dTriangulated.Points(convexHullIndices,:)';
  % Use geom3d trim extra points
  % result = convexHull(result);
end

function intersection_vertices = intersection(vertices1, vertices2)
  uniqueCols = @(array) unique(array', 'rows', 'stable')';
  vertices1 = uniqueCols(vertices1);
  vertices2 = uniqueCols(vertices2);
  poly1 = polyshape(vertices1(1,:), vertices1(2,:));
  poly2 = polyshape(vertices2(1,:), vertices2(2,:));
  poly_intersection = intersect(poly1, poly2);
  intersection_vertices = poly_intersection.Vertices';
end