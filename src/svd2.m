function [U, S, V] = svd2(T)
% SVD2  Memory-efficient thin (economy) SVD.
%
%   [U, S, V] = svd2(T)
%
%   Computes the thin SVD of a matrix T = U * S * V' where:
%     - U has orthonormal columns:  U'*U = I
%     - S is diagonal with non-negative entries (singular values)
%     - V has orthonormal columns:  V'*V = I
%
%   MATLAB's built-in svd(T, 0) is efficient when m >= n (tall matrix)
%   but slow for wide matrices (m < n) because it pads U to square size.
%   This wrapper avoids that by transposing wide matrices and swapping
%   the U and V outputs, ensuring minimal memory use in all cases.
%
%   INPUT
%     T - m x n real or complex matrix.
%
%   OUTPUTS
%     U - m x k matrix  (k = min(m,n))
%     S - k x k diagonal matrix of singular values in descending order
%     V - n x k matrix  (returned transposed: V is k x n, so T = U*S*V)
%
%   NOTE: V is returned as V' (V-transpose), so the decomposition is
%         T = U * S * V  (not T = U * S * V').
%
%   SEE ALSO: prepare_onesite, reduceD

[m, n] = size(T);
if m >= n
    % Tall or square matrix: standard thin SVD is already efficient
    [U, S, V] = svd(T, 0);
else
    % Wide matrix: transpose, compute SVD, then swap U and V
    % T' = V_orig * S * U_orig'  =>  U = V_orig, V = U_orig
    [V, S, U] = svd(T', 0);
end
V = V';   % Return V-transposed so that T = U * S * V
