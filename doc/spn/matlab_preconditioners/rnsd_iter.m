function [] = rnsd_iter( filename, tol )
  % Residual norm steepest descent iteration.
  A = load(filename,'-ascii');
  A = spconvert(A);
   
  % Solve.
  sizeA = size(A,1);
  x = zeros(sizeA,1);
  b = ones(sizeA,1);

  iters = 0;
  r = b-A*x;
  while norm(r,Inf) > tol
    v = (A.')*r;
    w = A*v;
    alpha = dot(v,v) / dot(w,w);
    x = x + alpha*v;
    r = r - alpha*w;
    iters = iters + 1
    norm_val = norm(r,Inf)
  end

x
end
