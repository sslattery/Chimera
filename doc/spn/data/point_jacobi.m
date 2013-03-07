function [] = point_jacobi( matname )
% Point-Jacobi preconditioned iteration matrix eigenvalues
filename = strcat(matname,'.mat');
A = load(filename,'-ascii');
A = spconvert(A); 
d = diag(A); 
d_inv = d.^(-1); 
M = diag(d_inv); 
I = speye(size(d_inv,1)); 
H = I-M*A;
opts.tol=1.0e-8;
eigs(H,1,'lm',opts)

end

