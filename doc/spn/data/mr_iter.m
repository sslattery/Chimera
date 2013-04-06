function [] = mr_iter( filename, tol )
  % MR iteration.
  A = load(filename,'-ascii');
  A = spconvert(A);
   
  % Solve.
  sizeA = size(A,1);
  x = zeros(sizeA,1);
  b = ones(sizeA,1);

  iters = 0;
  r = b-A*x;
  while norm(r,Inf) > tol
    r = b-A*x;
    Ar = A*r;
    alpha = dot(r,Ar) / dot(Ar,Ar);
    x = x + alpha*r;
    iters = iters + 1
    norm_val = norm(r,Inf)
  end

x
end