function [B, U, DB] = prepare_onesite(A, direction)
% PREPARE_ONESITE  Gauge-fix a single MPS tensor using SVD.
%
%   [B, U, DB] = prepare_onesite(A, direction)
%
%   Performs a QR-like decomposition (via thin SVD) on one MPS site tensor A
%   to bring it into left- or right-canonical form. This is the fundamental
%   step in MPS canonicalization.
%
%   An MPS tensor A has shape [D1, D2, d] where:
%     - D1 = left bond dimension
%     - D2 = right bond dimension
%     - d  = physical dimension
%
%   LEFT-CANONICAL ('lr'):
%     Reshapes A as (d*D1) x D2, computes A = B_matrix * (S*V'),
%     then reshapes B_matrix back to [D1, DB, d].
%     B satisfies B†B = I (left-canonical condition).
%     The matrix U = S*V' (shape [DB, D2]) is passed to the right neighbor.
%
%   RIGHT-CANONICAL ('rl'):
%     Reshapes A as D1 x (d*D2), computes A = (U*S) * B_matrix,
%     then reshapes B_matrix back to [DB, D2, d].
%     B satisfies BB† = I (right-canonical condition).
%     The matrix U = U*S (shape [D1, DB]) is passed to the left neighbor.
%
%   INPUTS
%     A         - MPS site tensor of shape [D1, D2, d].
%     direction - 'lr' for left-canonical, 'rl' for right-canonical.
%
%   OUTPUTS
%     B  - Isometric tensor in the specified canonical form.
%     U  - Singular value matrix to be absorbed into the neighboring site.
%     DB - New bond dimension after truncation (= rank of the reshaped A).
%
%   SEE ALSO: prepare, svd2, reduceD

[D1, D2, d] = size(A);

switch direction
    case 'lr'
        % Reshape to (d*D1) x D2: group physical+left into rows
        A = permute(A, [3, 1, 2]);         % [d, D1, D2]
        A = reshape(A, [d*D1, D2]);

        % Thin SVD: A = B_mat * S * V'
        [B, S, U] = svd2(A);
        DB = size(S, 1);                   % rank (new bond dimension)

        % Reshape left factor back to MPS tensor
        B = reshape(B, [d, D1, DB]);
        B = permute(B, [2, 3, 1]);         % [D1, DB, d]  (left-canonical)

        U = S * U;                         % [DB, D2]: absorbed into right site

    case 'rl'
        % Reshape to D1 x (d*D2): group physical+right into columns
        A = permute(A, [1, 3, 2]);         % [D1, d, D2]
        A = reshape(A, [D1, d*D2]);

        % Thin SVD: A = U * S * B_mat
        [U, S, B] = svd2(A);
        DB = size(S, 1);                   % rank (new bond dimension)

        % Reshape right factor back to MPS tensor
        B = reshape(B, [DB, d, D2]);
        B = permute(B, [1, 3, 2]);         % [DB, D2, d]  (right-canonical)

        U = U * S;                         % [D1, DB]: absorbed into left site
end
