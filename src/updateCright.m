function [Cright] = updateCright(Cright, B, X, A)
% UPDATECRIGHT  Extend the right environment by one site.
%
%   Cright = updateCright(Cright, B, X, A)
%
%   Updates the right boundary tensor (environment) by incorporating one
%   additional site moving right-to-left. This is the mirror operation of
%   updateCleft: instead of sweeping left-to-right, we sweep right-to-left.
%
%   The right environment Cright represents the partial contraction of the
%   network from the current site to site N:
%
%       Cright_{a, w, a'} = sum_{s,s'} A_{a, a_old, s}
%                                     * X_{w, w_old, s, s'}
%                                     * conj(B_{a', a'_old, s'})
%                                     * Cright_{a_old, w_old, a'_old}
%
%   After this call, Cright has shape [Dl_A, Dw_X_left, Dl_B] — the left
%   virtual indices of each layer, ready for use at the next site to the left.
%
%   INPUTS
%     Cright - Current right environment tensor of shape [Dr_A, Dw_X_r, Dr_B].
%              On the rightmost site (N+1) this is the scalar 1.
%     B      - Bra MPS tensor at the current site, shape [Dl_B, Dr_B, d].
%     X      - MPO tensor at the current site, shape [Dw_l, Dw_r, d, d].
%              Pass [] to use the identity operator.
%     A      - Ket MPS tensor at the current site, shape [Dl_A, Dr_A, d].
%
%   OUTPUT
%     Cright - Updated right environment, shape [Dl_A, Dw_l, Dl_B].
%
%   SEE ALSO: updateCleft, initCstorage

d = size(B, 3);
if isempty(X)
    % No operator: use identity MPO [1,1,d,d] = delta_{s,s'}
    X = reshape(eye(d), [1, 1, d, d]);
end

% Step 1: Contract ket A (index 2 = right bond) into Cright (index 3 = ket slot)
%   Result shape: [Dl_A, Dw_X_r*, Dr_B*, d_A]
Cright = contracttensors(A, 3, 2, Cright, 3, 3);

% Step 2: Contract MPO X (indices 2,4 = right bond and ket physical)
%   (index 4 of X ~ ket physical,  index 2 of X ~ MPO right bond in Cright)
%   Result shape: [Dw_l, d_bra, Dl_A, Dw_X_r*]  -->  rearranged
Cright = contracttensors(X, 4, [2, 4], Cright, 4, [4, 2]);

% Step 3: Contract bra conj(B) (indices 2,3 = right bond and bra physical)
%   Result shape: [Dw_l, Dl_B, Dl_A]  which serves as the new Cright
Cright = contracttensors(conj(B), 3, [2, 3], Cright, 4, [4, 2]);
