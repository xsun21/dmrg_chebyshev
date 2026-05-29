function [value] = HSnorm(mpoX, mpoY)
% Hilbert-Schmidt inner product of two Matrix Product Operators.
%
%   INPUTS
%     mpoX - Cell array {1,L} of MPO tensors with shape [Dl, Dr, d, d].
%     mpoY - Cell array {1,L} of MPO tensors with shape [Dl, Dr, d, d].
%
%   OUTPUT
%     value - Scalar


L = length(mpoX);
d = size(mpoX{1}, 3);   % physical dimension

% Reshape each site tensor from [Dl, Dr, d, d] to [Dl, Dr, d^2]
% to merge the two physical indices for contraction
for i = 1:L
    Xl = size(mpoX{i}, 1);   Xr = size(mpoX{i}, 2);
    Yl = size(mpoY{i}, 1);   Yr = size(mpoY{i}, 2);
    mpoX{i} = reshape(mpoX{i}, [Xl, Xr, d*d]);
    mpoY{i} = reshape(mpoY{i}, [Yl, Yr, d*d]);
end

% Contract site by site from left to right.
value = 1;   % scalar boundary condition at the left edge
for i = 1:L
    % Contract left bond of X* into current boundary (index 1 of X with index 1 of value)
    value = contracttensors(value, 2, 1, conj(mpoX{i}), 3, 1);
    value = permute(value, [2, 1, 3]);
    value = contracttensors(value, 3, [2, 3], mpoY{i}, 3, [1, 3]);
end

