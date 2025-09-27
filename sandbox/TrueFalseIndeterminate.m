classdef TrueFalseIndeterminate < int8
  % There are some basic tests in
  % `     run_TrueFalseIndeterminate_test_script.m

  enumeration
    true          ( 1);
    false         ( 0);
    indeterminate (-1);
    % true         ; % (true);
    % false        ; % (false);
    % indeterminate; % (false);
  end % End of enumeration members block

  methods(Access=public) % Define class methods.
    % ⋘──────── Override "logical" so that only "TrueFalseIndeterminate.true" treated as "truthy"  ────────⋙
    function tf = logical(this)
      tf = (this == TrueFalseIndeterminate.true);
    end % End of function

    % Overload "&" operator
    function result = and(left, right)
      if left == TrueFalseIndeterminate.true && right == TrueFalseIndeterminate.true 
        result = TrueFalseIndeterminate.true;
      elseif left == TrueFalseIndeterminate.false || right == TrueFalseIndeterminate.false
        result = TrueFalseIndeterminate.false;
      else
        result = TrueFalseIndeterminate.indeterminate;
      end
    end % End of function

    % Overload "|" operator
    function result = or(left, right)
      if left == TrueFalseIndeterminate.true || right == TrueFalseIndeterminate.true 
        result = TrueFalseIndeterminate.true;
      elseif left == TrueFalseIndeterminate.indeterminate || right == TrueFalseIndeterminate.indeterminate
        result = TrueFalseIndeterminate.indeterminate;
      else
        result = TrueFalseIndeterminate.false;
      end
    end % End of function

    % Overload "~" operator
    function result = not(this)
      if this == TrueFalseIndeterminate.true  
        result = TrueFalseIndeterminate.false;
      elseif this == TrueFalseIndeterminate.false  
        result = TrueFalseIndeterminate.true;
      elseif this == TrueFalseIndeterminate.indeterminate  
        result = TrueFalseIndeterminate.indeterminate;
      else
        error("TrueFalseIndeterminate.not(): Unexpected case.");
      end
    end % End of function

    function string_rep = char(this)
      % This method defines the string representation that is 
      % inserted when using "%s" in a string format, such as 
      % fprintf("%s", polyhedron) or sprintf("%s", polyhedron); 
      % It does not change the "disp" output. 
      arguments(Output)
        string_rep char; % Ensure output is cast to char, even if you create a string.
      end
      if this == TrueFalseIndeterminate.true  
        string_rep = sprintf("TrueFalseIndeterminate.true");
      elseif this == TrueFalseIndeterminate.false  
        string_rep = sprintf("TrueFalseIndeterminate.false");
      elseif this == TrueFalseIndeterminate.indeterminate  
        string_rep = sprintf("TrueFalseIndeterminate.indeterminate");
      else
        error("TrueFalseIndeterminate.char(): Unexpected case.");
      end
    end % End of function

    function disp(this)
      fprintf('%s', this);
    end % End of function

    % Overload "+" operator
    function result = plus(left, right)
      if class(left) == "string"
        result = left + char(right);
        return
      elseif class(right) == "string"
        result = char(left) + right;
        return
      end
      error("Only implemented for adding to strings.");
    end % End of function

    function out = myMethod(this)
      out = "Hello!";
    end % End of function

  end % End methods block

end % End of class