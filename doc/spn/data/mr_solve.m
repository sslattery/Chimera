function [] = mr_solve( filename, guess_scale, drop_tol, fill_level, mr_iters, tol, relax_param, max_iters )
  % MR preconditioned Richardson iteration.
  A = load(filename,'-ascii');
  A = spconvert(A);
  sizeA = size(A,1);

  % Build the inverse matrix.
  I = speye(sizeA);
  N = guess_scale*I;
  for j = 1:sizeA
    for i = 1:mr_iters
      r = I(:,j) - A*N(:,j);
      Ar = A*r;
      alpha = dot(r,Ar) / dot(Ar,Ar);
      N(:,j) = N(:,j) + alpha*r;
      N(:,j) = apply_drop( A, r, N(:,j), I(:,j), drop_tol, fill_level );
    end
  end

  % Solve.
  N = spconvert(N);
  H = I-A*N;
  spy(H);
  opts.tol=1.0e-8;
  eigs(H,1,'lm',opts)
  x = zeros(sizeA,1);
  b = ones(sizeA,1);

  iters = 0;
  r = b-A*N*x;
  while norm(r,Inf) > tol && (iters < max_iters)
    x = x + relax_param*r;
    r = b-A*N*x;
    iters = iters + 1
    norm(r,Inf)
  end

N*b
end

function [mj] = apply_drop( A, r, mj, ej, drop_tol, fill_level )
    % Chow and Sadd drop scheme
    n = size(mj,1);
    ATr = (A.')*r;
    deatr = dot(ej,ATr);
    rho = -2*deatr*mj + mj.^2;
    rho_sorted = sort(rho,'descend');
    for i = 1:n
	    if (rho(i,1) < rho_sorted(fill_level,1)) || (rho(i,1) < drop_tol)
	        mj(i,1) = 0.0;
	    end
    end
end
