function [U, S, V] = svd2(T)
% Memory-efficient SVD.
%
%   INPUT
%     T - m x n real or complex matrix.
%
%   OUTPUTS
%     U - m x k matrix 
%     S - k x k diagonal matrix of singular values in descending order
%     V - n x k matrix 
%

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
