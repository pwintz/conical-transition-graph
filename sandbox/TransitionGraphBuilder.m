classdef TransitionGraphBuilder < handle
  
  methods
    
    % Constructor
    function this = Transitioner(conic_abstractions)
      
      if nargin() == 0
        % Define construction of empty objects.;
      end
      
    end % End of constructor



  end % End of methods block

  methods(Abstract) % Define class methods.
    start_cone_ndxs = allStartConeIndices()
  end % End methods block

end % End of class