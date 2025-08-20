function str = num2strNear1(x)
  % Map a double "x" to a string with a fixed length of 10 with formatting that makes it clear the value relative to 1. 
  assert(isnumeric(x));
  assert(isreal(x));
  assert(isscalar(x));
  str_width = 11;
  if abs(x - 1) < 0.1
    if x > 1
      str = sprintf("1 + %.2g", x - 1);
    elseif x < 1
      str = sprintf("1 - %.2g", 1 - x);
    else
      str = sprintf("%" + str_width + "s", "1 (exact)");
    end
  else
    str = sprintf("%.2g", x);
  end
  str = sprintf("%" + str_width + "s", str);
  assert(strlength(str) == str_width, "length(str)=%d must equal str_width=%d for x=%g and str=%s", strlength(str), str_width, x, str)
end