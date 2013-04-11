function [] = ilut( filename, drop_tol )
% Block-Jacobi preconditioned iteration matrix eigenvalues
A = load(filename,'-ascii');
A = spconvert(A);
sizeA = size(A,1);

% Build the block preconditoner matrix
setup.type = 'ilutp';
setup.droptol = drop_tol;
[L,U] = ilu(A,setup);

% Compute the spectral radius
I = speye(sizeA);
H = I-(L^(-1))*A*(U^(-1));
opts.tol=1.0e-8;
eigs(H,1,'lm',opts)

end
