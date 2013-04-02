function [] = point_jacobi_solve( matname, tol )
% Point-Jacobi preconditioned iteration matrix eigenvalues
filename = strcat(matname,'.mat');
A = load(filename,'-ascii');
A = spconvert(A); 
d = diag(A); 
d_inv = d.^(-1); 
M = diag(d_inv); 
I = speye(size(d_inv,1)); 
H = I-M*A;

x = zeros(size(H,1))
b = ones(size(H,1))

iters = 0
while max(M*(b-A*x)) > tol
    x = H*x + M*b;
    iters = iters + 1
    max(M*(b-A*x))
end

x

end

