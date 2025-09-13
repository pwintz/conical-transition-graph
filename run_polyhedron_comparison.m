

A_ineq = [-1, 0; 0, -1; 1, 1];
b_ineq = [0; 0; 1];
check_lcon2vert_results(A_ineq, b_ineq);

A_ineq = [-1, 0; 0, -1; 1, 1; -1, -1];
b_ineq = [0; 0; 1; -1];
check_lcon2vert_results(A_ineq, b_ineq);


A_ineq = [-1, 0; 0, -1; 1, 1; -1, -1];
b_ineq = [0; 0; 1; -1];
A_eq = [-1, -1];
b_eq = [0];
check_lcon2vert_results(A_ineq, b_ineq, A_eq, b_eq);


function  check_lcon2vert_results(A_ineq, b_ineq, varargin)
  theirs_failed = false;
  mine_failed   = false;
  try
    verts_theirs = lcon2vert(A_ineq, b_ineq, varargin{:});
  catch exception
    theirs_failed = true;
  end

  try
    verts_mind   = polyhedron.lcon2vert(A_ineq, b_ineq, varargin{:});
  catch exception
    mine_failed = true;
  end
  
  assert(theirs_failed == mine_failed);
  if theirs_failed
    warning("Both lcon2vert's produced errors.");
  end

  pwintz.assertions.assertEqual(verts_theirs, verts_mind);
end