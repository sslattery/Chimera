function [] = rnsd_solve( filename, guess_scale, drop_tol, fill_level, mr_iters, tol, relax_param, max_iters )
  % Residual norm steepest descent preconditioned Richardson iteration.
  A = load(filename,'-ascii');
  A = spconvert(A);
  sizeA = size(A,1);

  % Build the inverse matrix.
  I = speye(sizeA);
  N = guess_scale*I;
  AT = A.';
  for j = 1:sizeA
    for i = 1:mr_iters
      r = I(:,j) - A*N(:,j);
      v = AT*r;
      w = A*v;
      alpha = dot(v,v) / dot(w,w);
      N(:,j) = N(:,j) + alpha*v;
      N(:,j) = apply_drop( v, N(:,j), I(:,j), drop_tol, fill_level );
    end
  end

  % Solve.
  N = spconvert(N);
  H = I-relax_param*A*N;
  spy(H,'.')
  opts.tol=1.0e-8;
  eigs(H,1,'lm',opts)
  x = zeros(sizeA,1);
  b = ones(sizeA,1);

  iters = 0;
  r = b-A*N*x;
  while norm(r,Inf) > tol && (iters < max_iters)
    x = x + relax_param*r;
    r = b-A*N*x;
    iters = iters + 1;
    norm(r,Inf);
  end
  iters
  N*x
end

function [mj] = apply_drop( v, mj, ej, drop_tol, fill_level )
  % Chow and Sadd drop scheme
  n = size(mj,1);
  rho = mj.^2 - 2*dot(ej,v)*mj;
  rho_sorted = sort(rho,'descend');
  for i = 1:n
    if (rho(i,1) < rho_sorted(fill_level,1)) || (rho(i,1) < drop_tol)
      mj(i,1) = 0.0;
    end
  end
end
