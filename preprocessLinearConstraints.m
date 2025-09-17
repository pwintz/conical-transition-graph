function [A_ineq, b_ineq, A_eq, b_eq] = preprocessLinearConstraints(A_ineq, b_ineq, A_eq, b_eq)
  % ` runtests Test_preprocessLinearConstraints.m

  arguments(Input)
    A_ineq;
    b_ineq;
    A_eq = [];
    b_eq = double.empty(0, 1);
  end % End of Input arguments block
  
  if isempty(A_eq)
    A_eq = double.empty(0, size(A_ineq, 2));
  end

  x_dim = size(A_ineq, 2);
  n_inequality_constraints = size(A_ineq, 1);
  n_equality_constraints   = size(A_eq,   1);
  pwintz.assertions.assertSize(A_eq,   [  n_equality_constraints, x_dim]);
  pwintz.assertions.assertSize(b_ineq, [n_inequality_constraints, 1]);
  pwintz.assertions.assertSize(b_eq,   [  n_equality_constraints, 1]);

  % Create normalized augmented matrices, where each row in A_ineq and A_eq have norm 1. 
  A_ineq_row_norms = pwintz.arrays.rowNorms(A_ineq);
  A_eq_row_norms   = pwintz.arrays.rowNorms(A_eq);
  A_ineq  = A_ineq ./ A_ineq_row_norms;
  b_ineq  = b_ineq ./ A_ineq_row_norms;
  A_eq  = A_eq     ./ A_eq_row_norms;
  b_eq  = b_eq     ./ A_eq_row_norms;

  constraints_dotted = A_ineq * A_ineq';

  tol = 1e-12; % The value of this tolerance is important. 
      parallel_constraints = constraints_dotted >=  1 - tol;
  antiparallel_constraints = constraints_dotted <= -1 + tol;

  is_ineq_const_to_delete        = false(n_inequality_constraints, 1);
  is_ineq_const_to_make_eq_const = false(n_inequality_constraints, 1);
  for i_parallel_constraint = find(parallel_constraints)'
    [m, n] = ind2sub(size(parallel_constraints), i_parallel_constraint);
    if m <= n
      % The matrix is symmetric, so we only need to check the upper or lower triangle.
      % The diagonal is not interesting because it contains the inner product of each row with itself.
      continue
    end
    % pwintz.strings.format("Constraint normals a_%d = %.2f and a_%d = %.2f are parallel.\nb_%d = %.2f and b%d = %.2f", m, A_ineq(m, :), n, A_ineq(n, :), m, b_ineq(m, :), n, b_ineq(n, :))
    if abs(b_ineq(m) - b_ineq(n) ) < tol
      % fprintf("Found redundant constraints. Will remove m = %d", m);
      is_ineq_const_to_delete(m) = true;
    end
  end

  for i_antiparallel_constraint = find(antiparallel_constraints)'
    [m, n] = ind2sub(size(antiparallel_constraints), i_antiparallel_constraint);
    if m <= n
      % The matrix is symmetric, so we only need to check the upper or lower triangle.
      % The diagonal is not interesting because it contains the inner product of each row with itself.
      continue
    end
    % pwintz.strings.format("Constraint normals a_%d = %.2f and a_%d = %.2f are antiparallel.\nb_%d = %.2f and b%d = %.2f", m, A_ineq(m, :), n, A_ineq(n, :), m, b_ineq(m, :), n, b_ineq(n, :))
    if abs(b_ineq(m) + b_ineq(n) ) < tol
      % fprintf("Found antiparallel constraints. Will remove m=%d and n=%d\n", m, n)
      is_ineq_const_to_delete(m) = true;
      is_ineq_const_to_delete(n) = true;
      % Add ineq constraint n to eq constraints (m would also work, but not both). 
      is_ineq_const_to_make_eq_const(n) = true; 
    end
  end

  % is_ineq_const_to_delete
  % is_ineq_const_to_make_eq_const
  % Add rows to A_eq, b_eq.
  % ! Must do this before deleting them from A_ineq, b_ineq!
  A_eq = [A_eq; A_ineq(is_ineq_const_to_make_eq_const, :)];
  b_eq = [b_eq; b_ineq(is_ineq_const_to_make_eq_const)];

  % Delete rows from A_ineq, b_ineq
  A_ineq = A_ineq(~is_ineq_const_to_delete, :);
  b_ineq = b_ineq(~is_ineq_const_to_delete);
  
  % if isempty(A_eq) 
  %   A_eq = double.empty(0, x_dim)
  % end
  % if isempty(A_ineq) 
  %   A_ineq = double.empty(0, x_dim)
  % end


% ! I don't thing the following lines are needed anymore.
%       Ab_ineq = [A_ineq, b_ineq];
%       Ab_eq   = [A_eq, b_eq];
% 
%       % Remove duplicate rows
%       Ab_ineq = pwintz.arrays.uniqueRows(Ab_ineq, tolerance=1e-9);
% 
%       % Find any rows that where the negation of the row is present in the array.
%       [~, duplicate_row_ndxs, ~, unique_row_ndxs] = pwintz.arrays.duplicatedRows([Ab_ineq; -Ab_ineq], tolerance=1e-9);
% 
% 
%       % all_indices = sort([duplicate_row_ndxs; unique_row_ndxs])
% 
%       % unique_row_ndxs    = intersect(1:size(A_ineq, 1), unique_row_ndxs)
%       % duplicate_row_ndxs = intersect(1:size(A_ineq, 1), duplicate_row_ndxs)
% 
%       ineq_rows_to_add_to_eq_constraints    = intersect(1:size(Ab_ineq, 1), duplicate_row_ndxs);
%       ineq_rows_to_keep_in_ineq_constraints = intersect(1:size(Ab_ineq, 1), unique_row_ndxs);
%       
%       % Add rows from Ab_ineq that create equality constraints.
%       Ab_eq = [Ab_eq; Ab_ineq(ineq_rows_to_add_to_eq_constraints, :)];
% 
      % % ╭────────────────────────────────────────────────────────────────────╮
      % % │             Remove duplicate equality constraint rows.             │
      % % ╰────────────────────────────────────────────────────────────────────╯
      % % TODO: This will not remove redundant rows if they have opposite signs. 
      % Ab_eq = pwintz.arrays.uniqueRows(Ab_eq, tolerance=1e-9);
      % 
      % % Remove any rows from Ab_ineq that we didnt add 
      % Ab_ineq = Ab_ineq(ineq_rows_to_keep_in_ineq_constraints, :);

      % [Ab_ineq_unique, src_ndxs_cell, out_rows_ndxs_cell] = pwintz.arrays.uniqueRows(C)

end