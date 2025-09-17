classdef FlowSpecification < handle
  
  properties(SetAccess = immutable) % Define instance constants.
    flow_map_matrix;
    flow_set_cone_ndxs;
  end % End of immutable properties block

  methods
    % Constructor
    function this = FlowSpecification(flow_set_cone_ndxs, flow_map_matrix)
      arguments(Input)
        flow_set_cone_ndxs      (1, :) {pwintz.validators.mustBeIndexVector};
        flow_map_matrix    (:, :) {pwintz.validators.mustBeSquare};
      end % End of Input arguments block
      
      this.flow_map_matrix = flow_map_matrix;
      this.flow_set_cone_ndxs = flow_set_cone_ndxs;
    end % End of constructor
  end % End of methods block
end % End of class