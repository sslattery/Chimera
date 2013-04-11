function [] = mr_iter( filename, tol )
% Minimal Residual (MR) iteration.
  A = load(filename,'-ascii');
  A = spconvert(A);
   
  % Solve.
  sizeA = size(A,1);
  x = zeros(sizeA,1);
  b = ones(sizeA,1);

  iters = 0;
  r = b-A*x;
  while norm(r,Inf) > tol
    Ar = A*r;
    alpha = dot(r,Ar) / dot(Ar,Ar);
    x = x + alpha*r;
    r = b-A*x;
    iters = iters + 1
    norm_val = norm(r,Inf)
  end

x
end
