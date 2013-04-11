function [] = richardson_solve( filename, relax_param, tol )
% Richardson iteration matrix eigenvalues
A = load(filename,'-ascii');
A = spconvert(A); 
I = speye(size(A,1)); 

x = zeros(size(A,1),1);
b = ones(size(A,1),1);

iters = 0;
r = b-A*x;
while norm(r,Inf) > tol
    x = x + relax_param*r;
    r = b-A*x;
    iters = iters + 1
    norm(r,Inf)
end

x

end
