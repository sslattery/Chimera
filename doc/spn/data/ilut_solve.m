function [] = ilut_solve( matname, drop_tol, tol )
% Block-Jacobi preconditioned iteration matrix eigenvalues
filename = strcat(matname,'.mat');
A = load(filename,'-ascii');
A = spconvert(A);
sizeA = size(A,1);

% Build the block preconditoner matrix
setup.type = 'ilutp';
setup.droptol = drop_tol;
[L,U] = ilu(A,setup);

% Compute the iteration matrix.
I = speye(sizeA);
M = L^(-1)
N = U^(-1)
H = I-M*A*N;

x = zeros(size(H,1))
b = ones(size(H,1))

iters = 0
while max(M*(b-A*x)) > tol
    x = H*x + M*b;
    iters = iters + 1
    max(M*(b-A*x))
end

N*x

end
