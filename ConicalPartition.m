classdef ConicalPartition < handle
  %%% ConicalPartition - Create a class for defining a partition of R^n into cones. 
  %%% The cones are defined usings a "watermelon slicing" method using equal angles for all cones.
  %%% To keep things simple, we are starting with  

  %%%
  
  properties(SetAccess = immutable)
    % Define instance constants.
    A % Matrix of "\dot x = A x".
    dimension  % The dimension of the system considered
    nSlices % Number of divisions per dimension  
    derivativeSpaceAngle % The angle swept by each cone in derivative space along the rotation in each dimension.
    
  end
  properties(SetAccess = immutable, GetAccess=private)
    derivativeSpaceSliceAngles % The angle of each slice in derivative space.
    derivativeVertices
    stateVertices
    stateVertexAngles
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
      % ⋘──── Compute Property Values Derived from Given Values ────⋙
      angles = linspace(0, 2*pi, this.nSlices + 1);
      this.derivativeSpaceSliceAngles = angles(1:end-1);
      this.derivativeVertices = nan([2, this.nSlices]);
      this.stateVertices = nan([2, this.nSlices]);
      % this.derivativeVertexAngles = nan([1, this.nSlices]);
      this.stateVertexAngles = nan([1, this.nSlices]);
      for i = 1:this.nSlices
        this.derivativeVertices(:,i) = this.getDerivativeSpaceVertex(i);
        this.stateVertices(:,i) = this.getStateSpaceVertex(i);
        this.stateVertexAngles(i) = atan2(this.stateVertices(2,i), this.stateVertices(1,i));
      end
      
      % ╭──────────────────────────────────────────╮
      % │             Check Properties             │
      % ╰──────────────────────────────────────────╯
      assert(size(this.A,1) == size(this.A,2), "The matrix A must be square")
      assert(cond(this.A) < 1e9, "The matrix A must be well-conditioned but had a condition number of %8.2g", cond(this.A))
    end

    function v = getStateSpaceVertex(this, varargin)
      % getVertex - Given the multi-index of the vertex, passed as the arguments, returns the vertex as a column vector.
      assert(numel(varargin) == 1, "we only support 2D systems now, so only one index should be given")
      v = this.getDerivativeSpaceVertex(varargin{1});
      v = inv(this.A)*v;
      v = v / norm(v);
    end

    function v = getDerivativeSpaceVertex(this, varargin)
      % getVertex - Given the multi-index of the vertex, passed as the arguments, returns the vertex as a column vector.
      assert(numel(varargin) == 1, "we only support 2D systems for now, so only one index should be given")
      % Compute the index of the vertex from the given argument, modulo the number of slices. Since MATLAB uses 1-based indexing, we need to subtract 1 before computing the modulo and then add it afterward.
      vertex_index = mod(varargin{1} - 1, this.nSlices) + 1;
      vertex_angle = this.derivativeSpaceSliceAngles(vertex_index);
      v = [cos(vertex_angle); sin(vertex_angle)];

    end

    function theta = getStateSpaceVertexAngle(this, varargin)
      % getVertex - Given the multi-index of the vertex, passed as the arguments, returns the vertex as a column vector.
      assert(numel(varargin) == 1, "we only support 2D systems now, so only one index should be given")
      vertex_index = mod(varargin{1} - 1, this.nSlices) + 1;
      theta = this.stateVertexAngles(vertex_index);
    end
    
    function C = getCone(this, varargin)
      % Given an index, get the state space cone for a given multi-index, given as a polyhedron (the cone truncated to the unit sphere), which is encoded as a row of vertices that defines the polyhedron.
      assert(numel(varargin) == 1, "we only support 2D systems now, so only one index should be given")
      vertex_index = mod(varargin{1} - 1, this.nSlices) + 1;
      v1 = this.getStateSpaceVertex(vertex_index);
      v2 = this.getStateSpaceVertex(vertex_index+1);
      origin = [0; 0];
      C = [v1, origin, v2];
    end

    function D = getDerivativeCone(this, varargin)
      % Given an index, get the state space cone for a given multi-index, given as a polyhedron (the cone truncated to the unit sphere), which is encoded as a row of vertices that defines the polyhedron.
      assert(numel(varargin) == 1, "we only support 2D systems now, so only one index should be given")
      vertex_index = mod(varargin{1} - 1, this.nSlices) + 1;
      v1 = this.getDerivativeSpaceVertex(vertex_index);
      v2 = this.getDerivativeSpaceVertex(vertex_index+1);
      origin = [0; 0];
      D = [v1, origin, v2];
    end
  end
end