classdef JumpSpecification < handle
  
  properties(SetAccess = immutable) % Define instance constants.
    % start_mode_ndx;
    % end_mode_ndx;
    jump_set_cone_ndxs;
    jump_map_matrix;
  end % End of immutable properties block

  methods
    % Constructor
    % function this = umpSpecification(start_mode_ndx, end_mode_ndx, jump_map_matrix, jump_set_cone_ndxs)
    function this = JumpSpecification(jump_set_cone_ndxs, jump_map_matrix)
      arguments(Input)
        % start_mode_ndx {pwintz.validators.mustBeIndexScalar};
        % end_mode_ndx {pwintz.validators.mustBeIndexScalar};
        jump_set_cone_ndxs      (:, :) {pwintz.validators.mustBeIndexVector};
        jump_map_matrix    (:, :) {pwintz.validators.mustBeSquare};
      end % End of Input arguments block
      
      % this.start_mode_ndx = start_mode_ndx;
      % this.end_mode_ndx = end_mode_ndx;
      this.jump_map_matrix = jump_map_matrix;
      this.jump_set_cone_ndxs = jump_set_cone_ndxs;
    end % End of constructor
  end % End of methods block
end % End of class