function [] = point_jacobi_solve( filename, relax_param, tol )
% Point-Jacobi preconditioned iteration matrix eigenvalues
A = load(filename,'-ascii');
A = spconvert(A); 
d = diag(A); 
d_inv = d.^(-1); 
M = diag(d_inv); 
I = speye(size(d_inv,1)); 

x = zeros(size(A,1),1);
b = ones(size(A,1),1);

iters = 0;
r = M*(b-A*x);
while norm(r,Inf) > tol
    x = x + relax_param*r;
    r = M*(b-A*x);
    iters = iters + 1
    norm(r,Inf)
end

x

end

