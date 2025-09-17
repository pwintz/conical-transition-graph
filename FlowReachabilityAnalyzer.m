classdef FlowReachabilityAnalyzer
  % ` runtests TestFlowReachabilityAnalyzer

  properties(SetAccess=immutable, GetAccess = {?matlab.unittest.TestCase}) % Define instance constants.
    conical_partition ConicalPartition;
    flow_set_cone_ndxs (1, :) {pwintz.validators.mustBeIndexVector};
    flow_map_matrix    (:, :) {pwintz.validators.mustBeSquare};
    can_flow_from_vertex_into_cone (:, :) logical;
    can_flow_from_cone_to_vertex   (:, :) logical;
    can_flow_from_vertex_to_vertex (:, :) logical;

    % The (i,j) entry of directly_reachable_sets_from_ray_through_cone contains 
    % the set that is reachable from vertex i by flowing through cone j.
    directly_reachable_sets_from_ray_through_cone      (:,:) cell; % ConvexPolyhedron; 
    restricted_reachable_sets_from_unit_sphere_in_cone (1,:) cell; % ConvexPolyhedron;
    directly_reachable_sets_from_vertex_to_vertex      (:,:) cell; % ConvexPolyhedron; 
  end % End of properties block

  methods
    function this = FlowReachabilityAnalyzer(conical_partition, flow_set_cone_ndxs, flow_map_matrix)
      arguments(Input)
        conical_partition ConicalPartition;
        flow_set_cone_ndxs (1, :) {pwintz.validators.mustBeIndexVector};
        flow_map_matrix (:, :)    {pwintz.validators.mustBeSquare};
      end % End of Input arguments block
      
      pwintz.assertions.assertAllAreMembers(flow_set_cone_ndxs, conical_partition.cone_indices);
      
      this.conical_partition  = conical_partition;
      this.flow_set_cone_ndxs = flow_set_cone_ndxs;
      this.flow_map_matrix    = flow_map_matrix;

      % ╭────────────────────────────────────────────────────────────────────────────────╮
      % │          Reachability from Unit Sphere within Each Cone ("Gain Sets")          │
      % ╰────────────────────────────────────────────────────────────────────────────────╯
      % Compute the region in each cone C_i that can be reached by a solution to \dot x = A_c * C_i that starts in unit sphere is always in C_i.
      restricted_reachable_sets_from_unit_sphere_in_cone = cell(conical_partition.n_cones, 1);
      for cone_ndx = this.flow_set_cone_ndxs
        restricted_reachable_sets_from_unit_sphere_in_cone{cone_ndx} = this.computeGainSetInCone(cone_ndx);
      end
      this.restricted_reachable_sets_from_unit_sphere_in_cone = restricted_reachable_sets_from_unit_sphere_in_cone;

      % ╭───────────────────────────────────────────────────────────────╮
      % │             Cone-Vertex Reachability Calculations             │
      % ╰───────────────────────────────────────────────────────────────╯
      % ⋘──────── Construct a list of cones that can be reached from each vertex ────────⋙
      % ⋘──────── Construct a list of vertices that can be reached from each cone ────────⋙
      can_flow_from_vertex_into_cone = false(conical_partition.n_vertices, conical_partition.n_cones);
      can_flow_from_cone_to_vertex   = false(conical_partition.n_cones, conical_partition.n_vertices);
      
      % Only check flows starting at rays, since linear systems cannont flow away from the origin.
      for ray_ndx = conical_partition.ray_indices
        ray      = conical_partition.getRay(ray_ndx);
        % ray_cone = ConvexPolyhedralCone.fromRays(ray);
        flow_dir = flow_map_matrix * ray; % Flow direction at "ray".

        % If the cone is not in the flow set, then it's impossible to flow into it.
        
        for adjacent_cone_ndx = this.adjacentFlowCones(ray_ndx) 
          adj_cone = conical_partition.getCone(adjacent_cone_ndx);
          
          % Test if flows.
          can_flow_into_cone = adj_cone.atXDoesVPointInward(ray, flow_dir);
          if can_flow_into_cone
            can_flow_from_vertex_into_cone(ray_ndx, adjacent_cone_ndx) = true;
            % [min_gain, max_gain] = this.gainsFromVertexToCone(ray_ndx, adjacent_cone_ndx);
          else
            can_flow_from_cone_to_vertex(adjacent_cone_ndx, ray_ndx) = true;
            % [min_gain, max_gain] = this.gainsFromConeToVertex(ray_ndx, adjacent_cone_ndx);
          end
        end % End of "adjacent_cone_ndx" for loop.
      end % End "ray_ndx" for loop
      this.can_flow_from_vertex_into_cone = can_flow_from_vertex_into_cone;
      this.can_flow_from_cone_to_vertex   = can_flow_from_cone_to_vertex;
      
      % ╭────────────────────────────────────────────────────────────────╮
      % │             Direction(s) of flows at each boundary             │
      % ╰────────────────────────────────────────────────────────────────╯
      % ! In 3D, we need to determine flow across 2D  boundaries instead of vertices, but in 2D the boundaries are vertices.
%       can_flow_from_bnd_into_cone = zeros(conical_partition.n_boundaries, conical_partition.n_cones);
%       can_flow_from_cone_into_bnd = zeros(conical_partition.n_cones, conical_partition.n_boundaries);
%       for i_bnd_ndx = conical_partition.boundary_indices
%         [bnd, adjacent_cone_ndxs, ray_indices] = conical_partition.getBoundary(i_bnd_ndx);
% 
%         can_flow_from_vertex_into_cone(i_ray_ndx, adjacent_cone_ndx) = true;
%         can_flow_from_cone_to_vertex(i_ray_ndx, adjacent_cone_ndx) = true;
% 
%         can_flow_from_bnd_into_cone(i_bnd_ndx)
%         can_flow_from_cone_into_bnd(i_bnd_ndx)
% 
%         rays = bnd.rays;
%         for ray = rays 
% 
%         end
%       end
      
      % can_flow_from_vertex_into_cone
      % can_flow_from_vertex_into_cone(conical_partition.origin_index, :)
      % assert(all(~can_flow_from_vertex_into_cone(conical_partition.origin_index, :)), "No cones can be reached from the origin.");
      

      
      % ╭────────────────────────────────────────────────────────────────────────────────────────────╮
      % │             Construct Flow Graph for Intra-cone Reachability between Vertices             │
      % ╰────────────────────────────────────────────────────────────────────────────────────────────╯
      directly_reachable_sets_from_ray_through_cone = cell(conical_partition.n_vertices, conical_partition.n_cones);
      directly_reachable_sets_from_vertex_to_vertex = cell(conical_partition.n_vertices, conical_partition.n_vertices);
      can_flow_from_vertex_to_vertex = false(conical_partition.n_vertices, conical_partition.n_vertices);
      % directly_reachable_sets_from_ray_through_cone = ConvexPolyhedron.arrayOfEmpty(conical_partition.n_vertices, conical_partition.n_cones);
      
      for ray_ndx = conical_partition.ray_indices
        % cone_ndxs = find(this.can_flow_from_vertex_into_cone(ray_ndx, :));
        % if isempty(cone_ndxs)
        %   continue
        % end
        % assert(isscalar(cone_ndx), "Is it possible that a flow only points in from one vertex? Guess we'll see! cone_ndx=%s", mat2str(cone_ndx))

        ray    = conical_partition.getRay(ray_ndx);
        
        for cone_ndx = this.conesDirectlyReachableFromVertex(ray_ndx)
          reach_set = this.computeReachableSet(cone_ndx, Polytope.fromPoint(ray));
          if isempty(reach_set)
            error('The reach set from from ray_ndx=%d in cone_ndx=%d is empty, despite the flow direction pointing from the ray into the cone: %s', ray_ndx, cone_ndx, reach_set);
          end
          if reach_set.isUnbounded()
            err = pwintz.Exception("FlowTransitionGainDigraph:UNBOUNDED_REACH_SET", ...
              "The reachable set from ray %d\n%D\nin cone %d\n%D\nis unbounded:\n%D",...  
              ray_ndx, ray, cone_ndx, conical_partition.getCone(cone_ndx), reach_set ...  
            );
            throw(err);
          end
          directly_reachable_sets_from_ray_through_cone{ray_ndx, cone_ndx} = reach_set;
          
          % ⋘──────── Get the reachable set vertices excluding the initial ray ────────⋙
          % initial_ray_ndx_in_reach_set_vertices = pwintz.arrays.findColumnIn(ray, reach_set_vertices, tolerance=1e-3, verbose=false);
          % reach_set_vertices = reach_set_vertices(:, setdiff(1:end, initial_ray_ndx_in_reach_set_vertices));
        
          % ⋘──────── Find which vertices are reachable and compute gains ────────⋙
          ray_norms = vecnorm(reach_set.vertices); % Used to compute gains.
          normalized_reach_set_vertices = nrv(reach_set.vertices);
          for bnd_vertex_ndx = setdiff(conical_partition.getVerticesAdjacentToCone(cone_ndx), ray_ndx)
            bnd_vertex = conical_partition.getVertex(bnd_vertex_ndx);
            
            col_ndxs = pwintz.arrays.findColumnIn(bnd_vertex, normalized_reach_set_vertices, tolerance=1e-5, verbose=false);
  
            if ~isempty(col_ndxs)
              can_flow_from_vertex_to_vertex(ray_ndx, bnd_vertex_ndx) = true;
              directly_reachable_sets_from_vertex_to_vertex{ray_ndx, bnd_vertex_ndx} = reach_set.vertices(:, col_ndxs);
            end
          end
        end
      end % End for loop
      this.directly_reachable_sets_from_ray_through_cone = directly_reachable_sets_from_ray_through_cone;
      this.directly_reachable_sets_from_vertex_to_vertex = directly_reachable_sets_from_vertex_to_vertex;
      this.can_flow_from_vertex_to_vertex = can_flow_from_vertex_to_vertex;
      assert(~isempty(directly_reachable_sets_from_ray_through_cone));
      assert(~isempty(directly_reachable_sets_from_vertex_to_vertex));
    end % End of constructor

    function adjacent_flow_set_cone_ndxs = adjacentFlowCones(this, ray_ndx) 
      arguments(Output)
        adjacent_flow_set_cone_ndxs (1, :);
      end % End of Input arguments block
      
      adjacent_flow_set_cone_ndxs = intersect(this.conical_partition.getConesAdjacentToRay(ray_ndx), this.flow_set_cone_ndxs);
    end

    function [min_gain, max_gain] = gainsFromVertexToCone(this, start_vertex_ndx, end_cone_ndx)
      if this.can_flow_from_vertex_into_cone(start_vertex_ndx, end_cone_ndx)
        % We use a  min_gain = max_gain = 1.0 when flowing from a vertex to a cone because the distance flowed to reach the cone is zero.  
        min_gain = 1.0;
        max_gain = 1.0;
      else 
        error("Cone not reachable.");
      end
    end

    function [min_gain, max_gain] = gainsFromConeToVertex(this, start_cone_ndx, end_vertex_ndx)
      if ~this.can_flow_from_cone_to_vertex(start_cone_ndx, end_vertex_ndx)
        pwintz.error(...
          "The vertex %d is not reachable from cone %d\nthis.can_flow_from_cone_to_vertex: %D\nverticesDirectlyReachableFromCone: %D", ...
          end_vertex_ndx, start_cone_ndx, this.can_flow_from_cone_to_vertex, this.verticesDirectlyReachableFromCone(start_cone_ndx)...
        )
      end

      adjacent_cone_ndxs = this.conical_partition.getConesAdjacentToVertex(end_vertex_ndx);
      pwintz.assertions.assertAllAreMembers(start_cone_ndx, adjacent_cone_ndxs);
      % 
      reachable_polytope = this.restricted_reachable_sets_from_unit_sphere_in_cone{start_cone_ndx};
      
      end_ray      = this.conical_partition.getRay(end_vertex_ndx);
      end_ray_cone = ConvexPolyhedralCone.fromRays(end_ray);
      is_reach_poly_verts_in_ray = end_ray_cone.containsPoints(reachable_polytope.vertices);
      
      % vert_ndxs_in_ray = pwintz.arrays.findColumnIn(end_ray, pwintz.arrays.normalizeColumns(reachable_polytope.vertices))
      reach_poly_verts_in_ray = reachable_polytope.vertices(:, is_reach_poly_verts_in_ray);
      gains = pwintz.arrays.columnNorms(reach_poly_verts_in_ray);
      min_gain = min(gains);
      max_gain = max(gains);

      % ! Debugging info
      % start_cone = this.conical_partition.getCone(start_cone_ndx);
      % reachable_polytope
      % end_ray      
      % end_ray_cone
      % is_reach_poly_verts_in_ray
      % reach_poly_verts_in_ray

      assert(~isempty(min_gain), "min_gain was empty");
      assert(~isempty(max_gain), "max_gain was empty");
      
    end
    
    function [min_gain, max_gain] = gainsFromVertexToVertex(this, start_vertex_ndx, end_vertex_ndx)
      reachable_set_vertices = this.directly_reachable_sets_from_vertex_to_vertex{start_vertex_ndx, end_vertex_ndx};
      if isempty(reachable_set_vertices)
        error("Vertex not reachable.");
      end
        
      % pwintz.strings.format("reachable_set vertices: %D", reachable_set_vertices);
      gains = pwintz.arrays.columnNorms(reachable_set_vertices);
      max_gain = max(gains);
      min_gain = min(gains);
      % Add transition from the non-origin "vertex_ndx" to the adjacent vertex "bnd_vertex_ndx" that can be reached by a single interval of flow.
      % this.addEdgeFromVertexToVertex(ray_ndx, bnd_vertex_ndx, min_gain, max_gain);
    end

    function reachable_cone_ndxs = conesDirectlyReachableFromVertex(this, vertex_ndx)
      reachable_cone_ndxs = find(this.can_flow_from_vertex_into_cone(vertex_ndx, :));
    end % End of function

    function reachable_vertex_ndxs = verticesDirectlyReachableFromCone(this, cone_ndx)
      reachable_vertex_ndxs = find(this.can_flow_from_cone_to_vertex(cone_ndx, :));
    end % End of function

    function reachable_vertex_ndxs = verticesDirectlyReachableFromVertex(this, vertex_ndx)
      reachable_vertex_ndxs = find(this.can_flow_from_vertex_to_vertex(vertex_ndx, :));
    end % End of function
  end % End of methods block.


  methods(Access = {?matlab.unittest.TestCase}) % Private block, accessible by unit tests.
    
    % ╭────────────────────────────────────────╮
    % │             Reachable Sets             │
    % ╰────────────────────────────────────────╯
    function [reachable_set_in_cone, reachable_set] = computeReachableSet(this, cone_ndx, P0)
      % Compute the cone that is reachable from P0 when flowing according to \dot x \in A*C_i, where C_i is given cone .
      % The "derivative cone" is represented using vertices on the unit sphere, so we scale it up so that we get a set the sufficiently represents "P0 + A*C_i", at least up to a large radius.
      arguments(Input)
        this;
        cone_ndx {pwintz.validators.mustBeIndexScalar};
        P0 Polytope;
      end % End of Input arguments block
      
      C_i  = this.conical_partition.getCone(cone_ndx);
      AC_i = this.flow_map_matrix * C_i;
    
      pwintz.assertions.assertAll(C_i.containsPoints(P0.vertices), "We have only implemented this function for the case where P0 is a subset of the cone.");
      
      % pwintz.strings.format("P0: %D", P0);
    
      % ⋘──────── Check for a degenerate case, w ────────⋙
      % If all of the vertices of P0 are in the boudnary of the cone, and AC_i points "out" of C_i, then the intersection is not full dimensional. Thus, it is hard to numerically find the intersection. 
      % In this case, the intersection will just be the portion of P0 that is in C_i. 
      do_any_point_inward = false;
      for v = this.flow_map_matrix * P0.vertices % For each direction
        for p = P0.vertices
          if C_i.atXDoesVPointInward(p, v)
            msg =pwintz.strings.format("At vertex p = %f, does the direction v = %f points into cone %d? %D", p, v, cone_ndx, C_i.atXDoesVPointInward(p, v));
            do_any_point_inward = true;
            break
          end
        end
      end
    

      reachable_set         = P0 + AC_i;

      if ~do_any_point_inward
        % If the cone does not point into the set, then it is numerically difficult to find the intersection, but it should just be the initial set.
        reachable_set_in_cone = P0;
        return;
      end
      
      % reachable_set_as_polytope = Polytope.fromConvexHull([P0.vertices + 100*AC_i.rays]);
      
      % % # Plot 
      % figure(1);
      % clf();
      % xlim(1.2*[-1, 1]);
      % ylim(1.2*[-1, 1]);
      % axis square;
      % hold on;
      % AC_i_as_polytope = AC_i.asPolytope();
      % C_i.plot("FaceColor", "red", "DisplayName", "C_i");
      % AC_i_as_polytope.plot("FaceColor", "blue", "DisplayName", "A*C_i");
      % P0.plot("FaceColor", "black", "DisplayName", "P0");
      % legend();
    
      try
        reachable_set_in_cone = intersection(C_i, reachable_set);
      catch err
        % If we can't compute the intersection, that probably means that flows point out of the cone, so the only points in the reachable set are the initial points.
        reachable_set_in_cone = P0;
      end
      
      % pwintz.strings.format("do_any_point_inward = %d", do_any_point_inward)
      % if do_any_point_inward
      %   % reachable_set_in_cone = intersection(reachable_set, C_i);
      % else
      %   reachable_set = P0;
      % end
    
      
      
    end
    
    
    function gain_set = computeGainSetInCone(this, cone_ndx)
      % For a given cone, compute the "gain set", which is the set within the cone that is reachable from the unit circle.
      sphere_overapproximation = this.conical_partition.getUnitSphereOverApproximationInCone(cone_ndx);
      [restricted_reachable_set, ~] = this.computeReachableSet(cone_ndx, sphere_overapproximation);
      assert(~isempty(restricted_reachable_set));
      gain_set = restricted_reachable_set;
    end
  end % End private methods block

end % End of class