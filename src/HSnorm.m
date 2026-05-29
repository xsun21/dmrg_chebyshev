function [value] = HSnorm(mpoX, mpoY)
% HSNORM  Hilbert-Schmidt inner product of two Matrix Product Operators.
%
%   value = HSnorm(mpoX, mpoY)
%
%   Computes the Hilbert-Schmidt (Frobenius) inner product between two
%   MPOs X and Y defined on the same lattice:
%
%       <X, Y>_HS = Tr[ X† Y ]
%
%   where the trace is over the full many-body Hilbert space. This is used
%   in ChebyshevH.m to estimate the truncation error after compressing the
%   Chebyshev MPO:
%
%       relative_error ≈ 1 - <T_trunc, T_full>_HS / <T_full, T_full>_HS
%
%   The computation is done by contracting the MPO network site by site,
%   maintaining a left-to-right boundary tensor of shape [Dr_X, Dr_Y].
%
%   INPUTS
%     mpoX - Cell array {1,L} of MPO tensors with shape [Dl, Dr, d, d].
%     mpoY - Cell array {1,L} of MPO tensors with shape [Dl, Dr, d, d].
%            Both MPOs must be defined on the same lattice with the same d.
%
%   OUTPUT
%     value - Scalar Tr[X† Y].
%
%   NOTE
%     When mpoX == mpoY this gives ||X||^2_HS (squared Frobenius norm).
%
%   SEE ALSO: ChebyshevH, contracttensors

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
% After site i, `value` is the partial contraction tensor of shape [Dr_X, Dr_Y].
% The index structure at each step is:
%   value_{a,b} = sum_{phys} conj(X_{...,a,phys}) * Y_{...,b,phys}
value = 1;   % scalar boundary condition at the left edge
for i = 1:L
    % Contract left bond of X* into current boundary (index 1 of X with index 1 of value)
    value = contracttensors(value, 2, 1, conj(mpoX{i}), 3, 1);
    % value now has shape: [Dr_value_left, Dr_X, phys]  ->  permute to [Dr_X, Dr_value_left, phys]
    value = permute(value, [2, 1, 3]);
    % Contract physical + left bond of Y into boundary: indices [2,3] of value with [1,3] of Y
    value = contracttensors(value, 3, [2, 3], mpoY{i}, 3, [1, 3]);
    % value now has shape [Dr_X, Dr_Y] — the new left boundary
end
% At the end, value is a scalar = Tr[X† Y]
