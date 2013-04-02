function [] = block_jacobi_solve( filename, spn, groups, tol )
% Block-Jacobi preconditioned iteration matrix eigenvalues
A = load(filename,'-ascii');
A = spconvert(A);
sizeA = size(A,1);

% Build the block preconditoner matrix
M = zeros(sizeA,sizeA);
block_size = groups*(spn+1)/2;
num_blocks = size(A,1)/block_size;
block_id = 0;
local_block = zeros(block_size,block_size);

for n = 1:num_blocks
    % Extract the block
    col_start = block_size*block_id;
    for i = 1:block_size
        row = col_start+i;
        for j = 1:block_size
            local_block(i,j) = A(row,col_start+j);
        end
    end
    
    % Invert the block
    local_block_inv = local_block^(-1);
    
    % Add it to the preconditioner
    col_start = block_size*block_id;   
    for i = 1:block_size
        row = col_start+i;
        for j = 1:block_size
            M(row,col_start+j) = local_block_inv(i,j);
        end
    end  
    
    block_id = block_id+1;
end

I = speye(sizeA);
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


