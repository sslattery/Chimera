function [] = richardson( filename, relax_param )
% Richardson iteration matrix eigenvalues
A = load(filename,'-ascii');
A = spconvert(A); 
I = speye(size(A,1)); 
H = I-relax_param*A;
opts.tol=1.0e-8;
eigs(H,1,'lm',opts)

end

