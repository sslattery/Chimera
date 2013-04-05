function [] = mr_solve( filename, guess_scale, drop_tol, fill_level mr_iters, tol )
  % MR preconditioned Richardson iteration.
  A = load(filename,'-ascii');
  A = spconvert(A);
  sizeA = size(A,1);


  % Build the inverse matrix.
  I = speye(sizeA);
  M = guess_scale*(speye(sizeA));
  for j = 1:sizeA
    M(:,j) = M*I(:,j);
    for i = 1:mr_iters
      r = I(:,j) - A*M(:,j);
      Ar = A*r;
      alpha = dot(r,Ar) / dot(Ar,Ar);
      M(:,j) = M(:,j) + alpha*r;
      apply_drop( (A.')*r, M(:,j), I(:,j), drop_tol, fill_level );
    end
  end


  % Solve.
  spy(I-A*M)
  x = zeros(sizeA,1);
  b = ones(sizeA,1);

  iters = 0;
  r = b-A*M*x;
  while norm(r,Inf) > tol
    x = x + relax_param*r;
    r = b-A*M*x;
    iters = iters + 1
    norm(r,Inf)
  end

M*x
end

function [] = apply_drop( ATr, mj, ej, drop_tol, fill_level )
    % Chow and Sadd drop scheme
    fill_counter = 0;
    n = size(mj,1);
    rho = -2*dot(ej,ATr)*mj + mj.^2;
    rho_sorted = sort(rho,'descend');
    for i = 1:n
	if rho[i] < rho[fill_level] || rho[i] < drop_tol
	   mj[i] = 0.0;
	end
    end
end
