function [Cright] = updateCright(Cright, B, X, A)
% Extend the right environment by one site.
%
%   INPUTS
%     Cright - Current right environment tensor of shape [Dr_A, Dw_X_r, Dr_B].
%     B      - Bra MPS tensor at the current site, shape [Dl_B, Dr_B, d].
%     X      - MPO tensor at the current site, shape [Dw_l, Dw_r, d, d].
%     A      - Ket MPS tensor at the current site, shape [Dl_A, Dr_A, d].
%
%   OUTPUT
%     Cright - Updated right environment, shape [Dl_A, Dw_l, Dl_B].
%

d = size(B, 3);
if isempty(X)
    % No operator: use identity MPO [1,1,d,d] = delta_{s,s'}
    X = reshape(eye(d), [1, 1, d, d]);
end

% Step 1: Contract ket A (index 2 = right bond) into Cright (index 3 = ket slot)
%   Result shape: [Dl_A, Dw_X_r*, Dr_B*, d_A]
Cright = contracttensors(A, 3, 2, Cright, 3, 3);

% Step 2: Contract MPO X (indices 2,4 = right bond and ket physical)
%   Result shape: [Dw_l, d_bra, Dl_A, Dw_X_r*] 
Cright = contracttensors(X, 4, [2, 4], Cright, 4, [4, 2]);

% Step 3: Contract bra conj(B) (indices 2,3 = right bond and bra physical)
%   Result shape: [Dw_l, Dl_B, Dl_A]  which serves as the new Cright
Cright = contracttensors(conj(B), 3, [2, 3], Cright, 4, [4, 2]);
