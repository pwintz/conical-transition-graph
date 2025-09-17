
cell_array = {1, 2, 3; 
              4, 5, 6};


% for i_row = 1:size(cell_array, 1)
%   element = cell_array(i_row, :);
%   element
% end % End of for block

cell_array_1d = {[1; 2], 2, 3, 4}
cell_array = cell_array_1d;

assert(isvector(cell_array), "Expected cell_array to be 1-dimensional. Instead its size was %s.", mat2str(size(cell_array)));
for ndx = 1:numel(cell_array)
  element = cell_array{ndx};

end
