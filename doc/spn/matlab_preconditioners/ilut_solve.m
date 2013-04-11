function [] = ilut_solve( filename, drop_tol, tol, relax_param )
% ILUT preconditioned Richardson iteration
A = load(filename,'-ascii');
A = spconvert(A);
sizeA = size(A,1);

% Build the LU decomposition
setup.type = 'ilutp';
setup.droptol = drop_tol;
[L,U] = ilu(A,setup);

% Compute the preconditoners
I = speye(sizeA);
M = L^(-1);
N = U^(-1);
spy(I-M*N*A)
x = zeros(sizeA,1);
b = ones(sizeA,1);

iters = 0;
r = M*b-M*A*N*x;
while norm(r,Inf) > tol
    x = x + relax_param*r;
    r = M*b-M*A*N*x;
    iters = iters + 1
    norm(r,Inf)
end

N*x

end
