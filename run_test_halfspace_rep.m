A_ineq = [1, 0; 
          1, 0;
          0, 1; 
          0, 1; 
          0, 1];
b_ineq = [ 3; 
           3; 
           1; 
           1; 
          -1];



% ⋘────────── Execute ─────────⋙
h = HalfspaceRepresentation(A_ineq, b_ineq);
h.A_ineq
h.b_ineq

dim = h.ambient_dimension;
for x = randn(dim, 1000)
  h_contains_x = h.containsPoint(x);
  Ax_leq_b = all(A_ineq * x <= b_ineq);
  pwintz.assertions.assertEqual(h_contains_x, Ax_leq_b, reason=sprintf("x=%s, Ax = %s, Ax - b = %s, Ax - b = %s (original).", mat2str(x), mat2str(h.A_ineq* x), mat2str(h.A_ineq* x - h.b_ineq), mat2str(A_ineq * x - b_ineq)));
end