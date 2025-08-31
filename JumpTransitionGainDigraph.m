classdef JumpTransitionGainDigraph < TransitionGainDigraph

  properties
    % Define instance variables.
    jump_set_ndxs;
    jump_set_cone_ndxs;
    jump_set_image_cone_ndxs;
    jump_map_matrix;
  end % End of properties block

  methods
    function this = JumpTransitionGainDigraph(conical_partition, jump_set_cone_ndxs, jump_set_image_cone_ndxs, jump_map_matrix)
      this@TransitionGainDigraph(conical_partition);
      this.jump_set_cone_ndxs = jump_set_cone_ndxs;
      this.jump_set_image_cone_ndxs = jump_set_image_cone_ndxs;
      this.jump_map_matrix    = jump_map_matrix;

      % ⋘──────── Alternative: If jump map is invertible ────────⋙
      assert(cond(this.jump_map_matrix) < 1e9, "Jump map must be invertible.");
      Ad_invs = inv(this.jump_map_matrix);
      for jump_set_image_cone_ndx = this.jump_set_image_cone_ndxs
        preimage_cone = Ad_invs * conical_partition.getCone(jump_set_image_cone_ndx);
        mean_of_preimage_cone = mean(preimage_cone.vertices, 2);
        preimage_ndx = conical_partition.getConesContainingPoint(mean_of_preimage_cone);
        origin_ndxs = pwintz.arrays.findColumnIn([0; 0], preimage_cone.vertices, tolerance=1e-6);
        preimage_nonzero_vertices = preimage_cone.vertices(:, setdiff(1:end, origin_ndxs));
        w_min = min(1./vecnorm(preimage_nonzero_vertices));
        w_max = max(1./vecnorm(preimage_nonzero_vertices));
        this.addEdgeFromConeToCone(preimage_ndx, jump_set_image_cone_ndx, w_min, w_max);
      end % End of for block

      % ⋘──────── Alternative for when jump map is not invertible ────────⋙
      
      %       for jump_cone_ndx = this.jump_set_cone_ndxs
      %         jump_cone_ndx
      %         jump_cone = conical_partition.getCone(jump_cone_ndx);
      %         jump_image_cone = this.jump_map_matrix * jump_cone;
      %         normalized_jump_image_vertices = nrv(jump_image_cone.vertices);
      %         jump_image_cone_normalized = ConvexPolyhedron.fromConvexHull(normalized_jump_image_vertices);
      %         nonzero_vertices = jump_image_cone_normalized.removeVertex([0; 0]).vertices;
      %         jump_image_angles = mod(pwintz.math.atan2(nonzero_vertices), 2*pi);
      %         if pwintz.math.angleDiffCCW(jump_image_angles, index=1) <= pi
      %           start_angle = jump_image_angles(1);
      %           end_angle = jump_image_angles(2);
      %         else 
      %           start_angle = jump_image_angles(2);
      %           end_angle = jump_image_angles(1);
      %         end
      %         [intersecting_cone_ndxs, arc] = conical_partition.getConesIntersectingArc(start_angle, end_angle);
      % 
      %         % ! This is a overly conservative estimate of the max and min norm of the cone mapped by 
      %         w_min = min(vecnorm(jump_image_cone.vertices));
      %         w_max = max(vecnorm(jump_image_cone.vertices));
      %         for jump_set_image_cone_ndx = intersecting_cone_ndxs
      %           jump_set_image_cone_ndx
      %           jump_transition_graph.addEdgeFromConeToCone(jump_cone_ndx, jump_set_image_cone_ndx, w_min, w_max);
      %         end
      %       end
    end
  end
end