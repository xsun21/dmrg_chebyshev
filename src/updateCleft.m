function [Cleft] = updateCleft(Cleft, B, X, A)
% UPDATECLEFT  Extend the left environment by one site.
%
%   Cleft = updateCleft(Cleft, B, X, A)
%
%   Updates the left boundary tensor (environment) by incorporating one
%   additional site moving left-to-right. The new Cleft represents the
%   partial contraction of the network from site 1 to the current site:
%
%       Cleft_{a', w', a} = sum_{s,s'} A_{a_old, a, s}
%                                    * X_{w_old, w', s, s'}
%                                    * conj(B_{a'_old, a', s'})
%                                    * Cleft_{a_old, w_old, a'_old}
%
%   The tensor index convention follows the MPS/MPO sandwich:
%
%       [bra B] ---Cleft--- [ket A]
%                  [MPO X]
%
%   After this call, Cleft has shape [Dr_A, Dw_X_right, Dr_B] — the right
%   virtual indices of each layer.
%
%   INPUTS
%     Cleft - Current left environment tensor of shape [Dl_A, Dw_X, Dl_B].
%             On the first site this is the scalar 1.
%     B     - Bra MPS tensor at the current site, shape [Dl_B, Dr_B, d].
%     X     - MPO tensor at the current site, shape [Dw_l, Dw_r, d, d].
%             Pass [] to use the identity operator (computes overlap <B|A>).
%     A     - Ket MPS tensor at the current site, shape [Dl_A, Dr_A, d].
%
%   OUTPUT
%     Cleft - Updated left environment, shape [Dr_A, Dw_r, Dr_B].
%
%   SEE ALSO: updateCright, initCstorage

d = size(B, 3);
if isempty(X)
    % No operator: use identity MPO [1,1,d,d] = delta_{s,s'}
    X = reshape(eye(d), [1, 1, d, d]);
end

% Step 1: Contract ket A (index 1 = left bond) into Cleft (index 3 = ket slot)
%   Result shape: [Dl_B*, Dw_X, Dr_A, d_A]   (via index 1 of A ~ index 3 of Cleft)
Cleft = contracttensors(A, 3, 1, Cleft, 3, 3);

% Step 2: Contract MPO X (indices 1,4 = left bond and ket physical) into Cleft
%   (index 4 of X ~ ket physical,  index 1 of X ~ index 2 of Cleft = MPO bond)
%   Result shape: [Dw_r, d_bra, Dl_B*, Dr_A]
Cleft = contracttensors(X, 4, [1, 4], Cleft, 4, [4, 2]);

% Step 3: Contract bra conj(B) (indices 1,3 = left bond and bra physical) into Cleft
%   (index 3 of B* ~ bra physical,  index 1 of B* ~ index 2 of Cleft = bra bond)
%   Result shape: [Dw_r, Dr_B, Dr_A]  ->  which is [Dr_A, Dw_r, Dr_B] after implicit permute
Cleft = contracttensors(conj(B), 3, [1, 3], Cleft, 4, [4, 2]);
