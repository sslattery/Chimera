function [] = richardson( matname )
% Richardson iteration matrix eigenvalues
filename = strcat(matname,'.mat');
A = load(filename,'-ascii');
A = spconvert(A); 
I = speye(size(A,1)); 
H = I-A;
opts.tol=1.0e-8;
eigs(H,1,'lm',opts)

end

