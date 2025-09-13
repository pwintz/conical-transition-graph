
for result = [TrueFalseIndeterminate.true, TrueFalseIndeterminate.false, TrueFalseIndeterminate.indeterminate]
result
  fprintf('=====\nValue: %s\n', result);
  fprintf('Negation: %s\n\n', ~result);
  if result
    disp(result + " is truthy");
  else
    disp(result  + " is not truthy");
  end
  if ~result
    disp("Negation of " + result + " is truthy");
  else
    disp("Negation of " + result  + " is not truthy");
  end
  % if 
end

TrueFalseIndederminate.true.myMethod()



% if TrueFalseIndeterminate.indeterminate
%   disp("This code is not called.");
% else
%   disp("This code is called.");
% end