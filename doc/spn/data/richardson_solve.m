function [] = richardson_solve( filename, relax_param, tol )
% Richardson iteration matrix eigenvalues
A = load(filename,'-ascii');
A = spconvert(A); 
I = speye(size(A,1)); 
H = I-relax_param*A;

x = zeros(size(H,1))
b = ones(size(H,1))

iters = 0
while max(b-A*x) > tol
    x = H*x + b;
    iters = iters + 1
    max(b-A*x)
end

x

end
