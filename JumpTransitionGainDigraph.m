classdef JumpTransitionGainDigraph < TransitionGainDigraph
  % ` runtests TestJumpTransitionGainDigraph

  properties(SetAccess=immutable)
    % Define instance variables.
    jump_set_cone_ndxs;
    jump_set_image_cone_ndxs;
    jump_map_matrix;
  end % End of properties block

  methods
    % function this = JumpTransitionGainDigraph(conical_partition, jump_set_cone_ndxs, jump_set_image_cone_ndxs, jump_map_matrix)
    function this = JumpTransitionGainDigraph(conical_partitions, jump_specifications)
      arguments(Input)
        conical_partitions  (:, 1) ConicalPartition;
        jump_specifications JumpSpecification {pwintz.validators.mustBeSquare};
        % jump_set_cone_ndxs (1, :) {pwintz.validators.mustBeIndexVector};
        % jump_map_matrices    double {pwintz.validators.mustBeSquare};
      end % End of Input arguments block
      
      conical_partitions_cell = mat2cell(conical_partitions, ones(numel(conical_partitions), 1));
      this@TransitionGainDigraph(conical_partitions_cell{:});

      jump_set_cone_ndxs       = cell(this.n_modes, this.n_modes);
      jump_set_image_cone_ndxs = cell(this.n_modes, this.n_modes);
      jump_map_matrix          = cell(this.n_modes, this.n_modes);


      for start_mode_ndx = 1:this.n_modes
        for end_mode_ndx = 1:this.n_modes
          jump_spec = jump_specifications(start_mode_ndx, end_mode_ndx);
          start_conical_partition = conical_partitions(start_mode_ndx);
          end_conical_partition   = conical_partitions(end_mode_ndx);

          pwintz.strings.format("start_conical_partition.n_cones: %d", start_conical_partition.n_cones)
          pwintz.strings.format("  end_conical_partition.n_cones: %d", end_conical_partition.n_cones)
          pwintz.strings.format("   jump_spec.jump_set_cone_ndxs: %d", jump_set_cone_ndxs)


          pwintz.assertions.assertAllAreMembers(jump_spec.jump_set_cone_ndxs, start_conical_partition.cone_indices);
          reach_analyzer = JumpReachabilityAnalyzer(start_conical_partition, end_conical_partition, jump_spec.jump_set_cone_ndxs, jump_spec.jump_map_matrix);

          jump_set_cone_ndxs{start_mode_ndx, end_mode_ndx} = jump_spec.jump_set_cone_ndxs;
          jump_set_image_cone_ndxs{start_mode_ndx, end_mode_ndx} = reach_analyzer.jumpImageSetConeIndices();
          jump_map_matrix{start_mode_ndx, end_mode_ndx} = jump_spec.jump_map_matrix;

          
          for start_cone_ndx = jump_spec.jump_set_cone_ndxs
            for end_cone_ndx = reach_analyzer.conesDirectlyReachableFromCone(start_cone_ndx)
              pwintz.strings.format("Can jump from Mode %d, Cone %d to Mode %d, Cone %d", start_mode_ndx, start_cone_ndx, end_mode_ndx, end_cone_ndx)

              [min_gain, max_gain] = reach_analyzer.gainsFromConeToCone(start_cone_ndx, end_cone_ndx);

              this.addEdgeFromConeToCone(start_mode_ndx, start_cone_ndx, end_mode_ndx, end_cone_ndx, min_gain, max_gain);
            end
          end
          % pwintz.assertions.assertAllAreMembers(jump_set_cone_ndxs, start_conical_partition.cone_indices);
        end % End of for block
      end % End of for block

      this.jump_set_cone_ndxs = jump_set_cone_ndxs;
      this.jump_set_image_cone_ndxs = jump_set_image_cone_ndxs;
      this.jump_map_matrix = jump_map_matrix;


      for start_mode_ndx = 1:this.n_modes
        for end_mode_ndx = 1:this.n_modes
          edges_between_modes = this.getEdgesFromConesToCones(start_mode_ndx, end_mode_ndx);
          if ~isempty(this.jump_set_cone_ndxs{start_mode_ndx})
            pwintz.assertions.assertNonempty(edges_between_modes);
          end
        end
      end
      return


%       % if ~isempty(jump_set_image_cone_ndxs)
%       %   pwintz.error("Passing jump_set_image_cone_ndxs to JumpTransitionGainDigraph is deprecated.")
%       % end
%       % this.jump_set_image_cone_ndxs = jump_set_image_cone_ndxs;
%       this.jump_map_matrix    = jump_map_matrix;
%       assert(cond(this.jump_map_matrix) < 1e9, "Jump map must be invertible.");
%       Ad_invs = inv(this.jump_map_matrix);
% 
%       is_jump_set_image_cone = false([1, conical_partition.n_cones]);
% 
%       % Find all the reachable cones for each jump cone.
%       for jump_set_cone_ndx = this.jump_set_cone_ndxs
%         jump_cone = conical_partition.getCone(jump_set_cone_ndx);
%         jump_angles = pwintz.math.atan2(jump_cone.rays);
%         jump_image_angles = mapAnglesByLinearMap(jump_map_matrix, jump_angles);
%         % jump_image = this.jump_map_matrix * jump_cone;
%         % jump_image_angles = pwintz.math.atan2(jump_image.rays);
% 
%         % The linear map can switch the order of the rays, as measured in the CCW direction, but the distance between them must be less than pi, since angle between the given rays is less than pi. 
%         if pwintz.math.angleDiffCCW(jump_image_angles, index=1) < pi
%           jump_image_start_angle = jump_image_angles(1);
%           jump_image_end_angle   = jump_image_angles(2);
%         else
%           jump_image_start_angle = jump_image_angles(2);
%           jump_image_end_angle   = jump_image_angles(1);
%         end
%           
%         % pwintz.math.angleDiffCCW(jump_image_angles_1)
%         pwintz.assertions.assertNumElements(jump_image_angles, 2);
%         pwintz.assertions.assertAllLessThan(pwintz.math.angleDiffCCW([jump_image_start_angle, jump_image_end_angle], index=1), pi, leftName="jump_image_angle difference (CCW)");
%         % jump_image_start_angle = jump_image_angles(1);
%         % jump_image_end_angle   = jump_image_angles(2);
%         image_cone_ndxs = conical_partition.getConesIntersectingArc(jump_image_start_angle, jump_image_end_angle);
%         is_jump_set_image_cone(image_cone_ndxs) = true;
%         for image_cone_ndx = image_cone_ndxs
%           image_cone_rays = conical_partition.getCone(image_cone_ndx).rays;
%           % To get the gains, we  
%           preimage_rays = Ad_invs * image_cone_rays;
%           w_min = min(1./vecnorm(preimage_rays));
%           w_max = max(1./vecnorm(preimage_rays));
%           this.addEdgeFromConeToCone(jump_set_cone_ndx, image_cone_ndx, w_min, w_max);
%         end 
% 
%         % mean_of_cone = mean(image_cone.rays, 2);
%         % image_ndx = conical_partition.getConesContainingPoint(mean_of_cone);
%         % % origin_ndxs = pwintz.arrays.findColumnIn([0; 0], preimage_cone.rays, tolerance=1e-6);
%         % % preimage_nonzero_rays = preimage_cone.rays(:, setdiff(1:end, origin_ndxs));
%         % w_min = min(vecnorm(image_cone.rays));
%         % w_max = max(vecnorm(image_cone.rays));
%         % this.addEdgeFromConeToCone(jump_set_cone_ndx, image_ndx, w_min, w_max);
%       end % End of for block
%       this.jump_set_image_cone_ndxs = find(is_jump_set_image_cone);

      % % ⋘──────── Alternative: If jump map is invertible ────────⋙
      % for jump_set_image_cone_ndx = this.jump_set_image_cone_ndxs
      %   preimage_cone = Ad_invs * conical_partition.getCone(jump_set_image_cone_ndx);
      %   mean_of_preimage_cone = mean(preimage_cone.rays, 2);
      %   preimage_ndx = conical_partition.getConesContainingPoint(mean_of_preimage_cone);
      %   origin_ndxs = pwintz.arrays.findColumnIn([0; 0], preimage_cone.rays, tolerance=1e-6);
      %   preimage_nonzero_rays = preimage_cone.rays(:, setdiff(1:end, origin_ndxs));
      %   w_min = min(1./vecnorm(preimage_nonzero_rays));
      %   w_max = max(1./vecnorm(preimage_nonzero_rays));
      %   this.addEdgeFromConeToCone(preimage_ndx, jump_set_image_cone_ndx, w_min, w_max);
      % end % End of for block

      % ⋘──────── Alternative for when jump map is not invertible ────────⋙
      
      %       for jump_cone_ndx = this.jump_set_cone_ndxs
      %         jump_cone_ndx
      %         jump_cone = conical_partition.getCone(jump_cone_ndx);
      %         jump_image_cone = this.jump_map_matrix * jump_cone;
      %         normalized_jump_image_rays = nrv(jump_image_cone.rays);
      %         jump_image_cone_normalized = ConvexPolyhedron.fromConvexHull(normalized_jump_image_rays);
      %         nonzero_rays = jump_image_cone_normalized.removeVertex([0; 0]).rays;
      %         jump_image_angles = mod(pwintz.math.atan2(nonzero_rays), 2*pi);
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
      %         w_min = min(vecnorm(jump_image_cone.rays));
      %         w_max = max(vecnorm(jump_image_cone.rays));
      %         for jump_set_image_cone_ndx = intersecting_cone_ndxs
      %           jump_set_image_cone_ndx
      %           jump_transition_graph.addEdgeFromConeToCone(jump_cone_ndx, jump_set_image_cone_ndx, w_min, w_max);
      %         end
      %       end
    end
  end
end