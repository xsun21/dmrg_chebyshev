function [Cleft] = updateCleft(Cleft, B, X, A)
% Extend the left environment by one site.
%
%   INPUTS
%     Cleft - Current left environment tensor of shape [Dl_A, Dw_X, Dl_B].
%     B     - Bra MPS tensor at the current site, shape [Dl_B, Dr_B, d].
%     X     - MPO tensor at the current site, shape [Dw_l, Dw_r, d, d].
%     A     - Ket MPS tensor at the current site, shape [Dl_A, Dr_A, d].
%
%   OUTPUT
%     Cleft - Updated left environment, shape [Dr_A, Dw_r, Dr_B].
%

d = size(B, 3);
if isempty(X)
    % No operator: use identity MPO [1,1,d,d] = delta_{s,s'}
    X = reshape(eye(d), [1, 1, d, d]);
end

% Step 1: Contract ket A (index 1 = left bond) into Cleft (index 3 = ket slot)
%   Result shape: [Dl_B*, Dw_X, Dr_A, d_A]   
Cleft = contracttensors(A, 3, 1, Cleft, 3, 3);

% Step 2: Contract MPO X (indices 1,4 = left bond and ket physical) into Cleft
%   Result shape: [Dw_r, d_bra, Dl_B*, Dr_A]
Cleft = contracttensors(X, 4, [1, 4], Cleft, 4, [4, 2]);

% Step 3: Contract bra conj(B) (indices 1,3 = left bond and bra physical) into Cleft
%   Result shape: [Dw_r, Dr_B, Dr_A]  
Cleft = contracttensors(conj(B), 3, [1, 3], Cleft, 4, [4, 2]);
